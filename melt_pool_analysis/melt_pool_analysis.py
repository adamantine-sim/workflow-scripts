import os
import glob
import re

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pyvista


def melt_pool_analysis(
    path_to_adamantine_files,
    adamantine_filename,
    output_directory,
    scalar_name="temperature",
    isovalue=1670.0,
    show_plot=True,
):
    """
    Perform melt pool analysis on Adamantine .pvtu/.vtu outputs.

    Parameters
    ----------
    path_to_adamantine_files : str
        Directory containing solution.<iter>.pvtu and solution.<iter>.<rank>.vtu files.
    adamantine_filename : str
        Prefix of the Adamantine output files, e.g. "solution".
    output_directory : str
        Directory where plots/results should be written.
    scalar_name : str, optional
        Active scalar field to use for clipping. Default is "temperature".
    isovalue : float, optional
        Temperature isovalue used to define the melt pool. Default is 1670.0.
    show_plot : bool, optional
        If True, display the plot interactively. Default is True.

    Returns
    -------
    dict
        Dictionary containing per-step arrays and average melt pool dimensions.
    """
    print("Performing melt pool analysis")

    data_dir = os.path.abspath(path_to_adamantine_files)
    output_dir = os.path.abspath(output_directory)
    os.makedirs(output_dir, exist_ok=True)

    time_series_plot_filename = os.path.join(output_dir, "melt_pool_ts_plot.png")

    # ----------------------------------------------------------
    # Plot styling
    # ----------------------------------------------------------
    SMALL_SIZE = 12
    MEDIUM_SIZE = 16
    BIGGER_SIZE = 18

    plt.rc("font", size=SMALL_SIZE)
    plt.rc("axes", titlesize=BIGGER_SIZE)
    plt.rc("axes", labelsize=MEDIUM_SIZE)
    plt.rc("xtick", labelsize=SMALL_SIZE)
    plt.rc("ytick", labelsize=SMALL_SIZE)
    plt.rc("legend", fontsize=SMALL_SIZE)
    plt.rc("figure", titlesize=BIGGER_SIZE)

    # ----------------------------------------------------------
    # Helpers
    # ----------------------------------------------------------
    def get_iteration_numbers():
        """
        Find iteration numbers from .pvtu files only.
        Example match: solution.123.pvtu -> iteration 123
        """
        pattern = os.path.join(data_dir, f"{adamantine_filename}.*.pvtu")
        files = glob.glob(pattern)

        iteration_numbers = []
        regex = re.compile(rf"^{re.escape(adamantine_filename)}\.(\d+)\.pvtu$")

        for f in files:
            base = os.path.basename(f)
            match = regex.match(base)
            if match:
                iteration_numbers.append(int(match.group(1)))

        return sorted(set(iteration_numbers))

    def extract_time_from_vtu(vtu_path):
        """
        Try to extract simulation time from a .vtu file.
        First tries to find a numeric value between XML tags on a line containing TimeValue.
        Falls back to the first numeric XML-tag value found after scanning the file.
        """
        if not os.path.isfile(vtu_path):
            return np.nan

        with open(vtu_path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()

        # Prefer lines that mention time explicitly
        for line in lines:
            if "TimeValue" in line or "time" in line.lower():
                m = re.search(r">([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)<", line)
                if m:
                    return float(m.group(1))

        # Fallback: preserve behavior similar to the older script
        if len(lines) > 8:
            m = re.search(r">([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)<", lines[8])
            if m:
                return float(m.group(1))

        return np.nan

    # ----------------------------------------------------------
    # Discover timesteps
    # ----------------------------------------------------------
    iteration_numbers = get_iteration_numbers()
    if not iteration_numbers:
        raise FileNotFoundError(
            f"No .pvtu files found for prefix '{adamantine_filename}' in '{data_dir}'. "
            f"Expected files like '{adamantine_filename}.<iteration>.pvtu'."
        )

    # ----------------------------------------------------------
    # Process each timestep
    # ----------------------------------------------------------
    cycle_index = []
    depth = []
    length = []
    width = []
    sim_time = []

    # No need to keep a persistent plotter around for this analysis
    for iteration_number in iteration_numbers:
        pvtu_path = os.path.join(data_dir, f"{adamantine_filename}.{iteration_number}.pvtu")
        vtu_rank0_path = os.path.join(data_dir, f"{adamantine_filename}.{iteration_number}.0.vtu")

        if not os.path.isfile(pvtu_path):
            print(f"Warning: missing {pvtu_path}; skipping")
            continue

        dataset = pyvista.read(pvtu_path)

        # Handle scalar name more safely
        available_scalars = list(dataset.array_names)
        chosen_scalar = None

        if scalar_name in available_scalars:
            chosen_scalar = scalar_name
        else:
            # common fallback if capitalization differs
            for name in available_scalars:
                if name.lower() == scalar_name.lower():
                    chosen_scalar = name
                    break

        if chosen_scalar is None:
            raise KeyError(
                f"Scalar '{scalar_name}' not found in {pvtu_path}. "
                f"Available scalars: {available_scalars}"
            )

        dataset.set_active_scalars(chosen_scalar)
        clipped_volume = dataset.clip_scalar(scalars=chosen_scalar, value=isovalue, invert=False)

        melt_pool_depth = 0.0
        width_val = 0.0
        length_val = 0.0

        n_points = len(clipped_volume.points)

        if n_points > 0:
            bounds = clipped_volume.bounds
            melt_pool_depth = abs(bounds[5] - bounds[4])

            # Project 3D points to XY for ellipse fit
            points = clipped_volume.points
            points_2d = points[:, :2].astype(np.float32).reshape(-1, 1, 2)

            # OpenCV fitEllipse requires at least 5 points
            if len(points_2d) >= 5:
                ellipse = cv2.fitEllipse(points_2d)
                (_, _), (axis_a, axis_b), _ = ellipse
                length_val = float(max(axis_a, axis_b))
                width_val = float(min(axis_a, axis_b))
            else:
                # Too few points for ellipse fitting; keep depth, zero out width/length
                width_val = 0.0
                length_val = 0.0

        cycle_index.append(iteration_number)
        depth.append(float(melt_pool_depth))
        width.append(float(width_val))
        length.append(float(length_val))
        sim_time.append(float(extract_time_from_vtu(vtu_rank0_path)))

    if not cycle_index:
        raise RuntimeError("No valid timesteps were processed.")

    # ----------------------------------------------------------
    # Plot and save
    # ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 7.5))
    ax.plot(sim_time, depth, label="depth")
    ax.plot(sim_time, width, label="width")
    ax.plot(sim_time, length, label="length")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Dimensions (m)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(time_series_plot_filename, dpi=300)

    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    # ----------------------------------------------------------
    # Summary stats
    # ----------------------------------------------------------
    depth_arr = np.array(depth, dtype=float)
    width_arr = np.array(width, dtype=float)
    length_arr = np.array(length, dtype=float)

    nonzero_depth = depth_arr[depth_arr > 0]
    nonzero_width = width_arr[width_arr > 0]
    nonzero_length = length_arr[length_arr > 0]

    avg_depth = float(np.mean(nonzero_depth)) if len(nonzero_depth) else 0.0
    avg_width = float(np.mean(nonzero_width)) if len(nonzero_width) else 0.0
    avg_length = float(np.mean(nonzero_length)) if len(nonzero_length) else 0.0

    print(f"Melt-pool plot: {time_series_plot_filename}")
    print(f"Avg meltpool width: {avg_width}")
    print(f"Avg meltpool depth: {avg_depth}")
    print(f"Avg meltpool length: {avg_length}")

    return {
        "plot_filename": time_series_plot_filename,
        "iteration": cycle_index,
        "time": sim_time,
        "depth": depth,
        "width": width,
        "length": length,
        "avg_width": avg_width,
        "avg_depth": avg_depth,
        "avg_length": avg_length,
    }
