import numpy as np
import pyvista
import os
import glob
import matplotlib.pyplot as plt
import sys
from scipy.spatial import KDTree
import argparse


# ----------------------------------------------------------
# Function to get the simulation iteration count (should be generalized across workflow-scripts)
# ----------------------------------------------------------

def get_iteration_count(path_prefix):
    filename_pattern = path_prefix + '*.pvtu'
    files = glob.glob(filename_pattern)

    iteration_numbers = []
    for f in files:
        temp = f.replace('.pvtu', '')
        post_id_index = temp.rfind('.')
        iteration_string = temp[post_id_index+1:]
        iteration_count = int(iteration_string)
        iteration_numbers.append(iteration_count)

    iteration_numbers.sort()

    return iteration_numbers

def check_indexing(directory, adamantine_filename, iteration_numbers):
    # This function is not currently used since the indexing is not consistent
    print("Checking indexing...")
    
    filename = directory + adamantine_filename + '.' + str(iteration_numbers[0]) + '.pvtu'
    mesh_0 = pyvista.read(filename)

    consistent_indexing = False

    for iteration_number in iteration_numbers:
        filename = directory + adamantine_filename + '.' + str(iteration_number) + '.pvtu'
        mesh_i = pyvista.read(filename)

        consistent_indexing = mesh_0.number_of_points == mesh_i.number_of_points and np.allclose(mesh_0.points, mesh_i.points)

        if not consistent_indexing:
            print("ERROR: Indexing is not consistent between time steps! Iteration number: ", iteration_number)
            sys.exit()
        
    return consistent_indexing

def get_time_series(directory, adamantine_filename, field_name, line_plots):
    # NOTE: About 1 out of every 10 points has one or more spurious drops to zero after being activated. I'm not sure what
    # that is caused by. Maybe this is due to points that lie on the activated/deactivated interface?

    # Would it be faster to sort each time series by location instead of finding these maps every time?

    print("Getting the time series...")

    iteration_numbers = get_iteration_count(os.path.join(directory, adamantine_filename))

    # Load first dataset (reference mesh)
    dataset_0 = pyvista.read(f"{directory}/{adamantine_filename}.{iteration_numbers[0]}.pvtu")
    n_points = dataset_0.number_of_points

    # Allocate the time series array
    val = np.zeros((len(iteration_numbers), n_points))

    # Preallocate tree query result (optional, for speed)
    distances = np.empty(n_points)
    nearest_indices = np.empty(n_points, dtype=int)


    ref_points = dataset_0.points

    for idx, iteration_number in enumerate(iteration_numbers):
        file_to_plot = f"{directory}/{adamantine_filename}.{iteration_number}.pvtu"
        dataset = pyvista.read(file_to_plot)

        # Build KDTree of current timestep points
        tree = KDTree(dataset.points)

        # Match reference points to nearest in current dataset
        distances, nearest_indices = tree.query(ref_points, workers=-1)

        if np.sum(distances) > 0.01:
            print(f'WARNING: timestep {iteration_number} has large total match distance: {np.sum(distances):.4f}')

        # Vectorized temperature lookup
        val[idx] = dataset.point_data[field_name][nearest_indices]

    # Optional plotting for to see the single-point temperature history
    if (line_plots):
        for i in range(min(5, val.shape[1])):  # plot first 5 points
            plt.figure()
            plt.plot(iteration_numbers, val[:, i], 'b.-')
            plt.title(f"Point {i} Temperature Over Time")
            plt.xlabel("Iteration")
            plt.ylabel("Temperature")
            plt.grid(True)
            plt.show()

    return val

def sort_points_by_xyz(points):
    # Use structured array to sort lexicographically by x, then y, then z
    return np.lexsort((points[:, 2], points[:, 1], points[:, 0]))

def get_time_series_sort(directory, adamantine_filename, field_name, line_plots):
    # NOTE: About 1 out of every 10 points has one or more spurious drops to zero after being activated. I'm not sure what
    # that is caused by. Maybe this is due to points that lie on the activated/deactivated interface?

    print("Getting the time series...")

    iteration_numbers = get_iteration_count(os.path.join(directory, adamantine_filename))

    # Load reference
    ref_dataset = pyvista.read(f"{directory}/{adamantine_filename}.{iteration_numbers[0]}.pvtu")
    ref_points = ref_dataset.points
    sort_idx = sort_points_by_xyz(ref_points)
    sorted_ref_points = ref_points[sort_idx]
    n_points = len(sorted_ref_points)

    # Initialize result
    val = np.empty((len(iteration_numbers), n_points))

    for idx, iteration_number in enumerate(iteration_numbers):
        file_path = f"{directory}/{adamantine_filename}.{iteration_number}.pvtu"
        dataset = pyvista.read(file_path)

        points = dataset.points
        if points.shape != ref_points.shape:
            raise ValueError(f"Point count mismatch at timestep {iteration_number}")

        # Sort this timestep's points
        current_sort_idx = sort_points_by_xyz(points)
        sorted_points = points[current_sort_idx]

        # Compare sorted points to reference
        if not np.allclose(sorted_points, sorted_ref_points, atol=1e-10):
            print(f"WARNING: timestep {iteration_number} point mismatch after sorting!")

        # Reorder data using sort index
        temperature = dataset.point_data[field_name][current_sort_idx]
        val[idx] = temperature

    # Optional plotting for to see the single-point temperature history
    if (line_plots):
        for i in range(min(5, val.shape[1])):  # plot first 5 points
            plt.figure()
            plt.plot(iteration_numbers, val[:, i], 'b.-')
            plt.title(f"Point {i} Temperature Over Time")
            plt.xlabel("Iteration")
            plt.ylabel("Temperature")
            plt.grid(True)
            plt.show()

    # TODO: Do I need to try to get these back into the "right" order for plotting?

    return val


def objective_17_4_PH(time_series):
    # NOTE: I have a special check in here to ignore cases where the temperature dips down below 5K to deal with artifacts
    # in the data. I'm still not sure why those occur.

    K_to_C = 273.15

    T_solidus = 1713.0 # Solidus temperature for 17-4PH
    T_Ms = 150.0 + K_to_C # Martensite start temperature in K
    T_ppt = 480.0 + K_to_C # Optimal temperature for precipitation
    T_band_ppt = 50.0 # Width of the band around T_ppt that "count" for the objective function, this is the "diameter" of the band. The band has a hard upper limit at the reaustenization temperature
    T_A = 550.0 + K_to_C # Re-austenization temperature

    per_point_scores = np.zeros(time_series.shape[1])
    
    # Step through time
    for n in range(0, time_series.shape[1]):
        single_point_series = time_series[:,n]

        has_melted = False
        is_martensite = False

        for temperature in single_point_series:
            if temperature >= T_solidus:
                    has_melted = True
                    per_point_scores[n] = 0
            elif has_melted and (not is_martensite):
                if temperature > 5.0 and temperature < T_Ms:
                    is_martensite = True
            elif is_martensite:
                if temperature > T_A:
                    is_martensite = False
                elif np.abs(temperature-T_ppt) < T_band_ppt/2.0:
                    per_point_scores[n] = per_point_scores[n] + 1

    return per_point_scores


def plot_score_on_mesh(directory, adamantine_filename, per_point_scores):
    iteration_numbers = get_iteration_count(directory + adamantine_filename)
    filename = directory + adamantine_filename + '.' + str(iteration_numbers[0]) + '.pvtu'
    dataset = pyvista.read(filename)

    dataset.point_data['score'] = per_point_scores
    pl = pyvista.Plotter()
    dataset.set_active_scalars("score")
    clipped = dataset.clip('y')
    pl.add_mesh(clipped, show_edges=True)
    pl.show()


def analysis(directory, adamantine_filename, field_name, line_plots, volume_plot):
    #time_series = get_time_series(directory, adamantine_filename, field_name, line_plots)
    time_series = get_time_series_sort(directory, adamantine_filename, field_name, line_plots)

    per_point_scores = objective_17_4_PH(time_series)

    if volume_plot:
        plot_score_on_mesh(directory, adamantine_filename, per_point_scores)

    # I may also want to promote low variance between the score at different locations
    score = np.sum(per_point_scores)

    print("Score:", score)

    return score

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Analysis routine for 17-4PH",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--data-directory", help="location of the VTK files to load")
    parser.add_argument("-f", "--base-filename", help="base filename for the VTK files")
    parser.add_argument("-v", "--variable", default="temperature", help='variable name in the VTK files')
    parser.add_argument("--line-plots", default=False, help='whether to plot some single-point temperature histories')
    parser.add_argument("--volume-plot", default=False, help='whether to plot a 3D cross-section of the per-point score')
    args = parser.parse_args()

    directory = args.data_directory #'/Users/71d/Documents/workspace/DockerWorkspace_adamantine/power_rook_2x/'
    adamantine_filename = args.base_filename #'solution_dwell'
    field_name = args.variable #"temperature"
    line_plots = args.line_plots
    volume_plot = args.volume_plot

    analysis(directory, adamantine_filename, field_name, line_plots, volume_plot)