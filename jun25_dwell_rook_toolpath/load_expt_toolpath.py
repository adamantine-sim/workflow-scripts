import csv
import copy
import sys
import os
from pathlib import Path
from typing import List

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from toolpath_writer import create_toolpath
from toolpath_writer import get_toolpath_info
from toolpath_writer import generate_print_plan_file
from toolpath_writer import write_event_series


def split_csv(input_file, output_dir):
    output_dir = Path(output_dir)
    print_dir = output_dir / "print_layers"
    reheat_dir = output_dir / "reheat_layers"

    print_dir.mkdir(parents=True, exist_ok=True)
    reheat_dir.mkdir(parents=True, exist_ok=True)

    epsilon = 1.0e-10

    # Counters for file naming (1-based indexing)
    print_idx = 1
    reheat_idx = 1

    # Buffers
    current_rows = []       # rows for the active block
    pending_zeros = []      # rows of zeros before the next block
    current_type = None     # 'print' or 'reheat'

    def compress_zeros(rows):
        """Compress stretches of zeros to keep only 0th, 1st, and last entry"""
        if not rows:
            return rows
        result = []
        zero_run = []
        for row in rows:
            value = float(row[4])
            if value == 0.0:
                zero_run.append(row)
            else:
                # end of a zero run
                if zero_run:
                    if len(zero_run) <= 2:
                        result.extend(zero_run)
                    else:
                        result.extend([zero_run[0], zero_run[1], zero_run[-1]])
                    zero_run = []
                result.append(row)
        # handle a run at the very end
        if zero_run:
            if len(zero_run) <= 2:
                result.extend(zero_run)
            else:
                result.extend([zero_run[0], zero_run[1], zero_run[-1]])
        return result

    def flush_rows(rows, out_type, idx):
        """Write accumulated rows to file"""
        if not rows:
            return
        rows = compress_zeros(rows)
        if out_type == "print":
            out_file = print_dir / f"layer_{idx}_scan_path.txt"
        else:
            out_file = reheat_dir / f"layer_{idx}_scan_path.txt"

        with open(out_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(rows)

    with open(input_file, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 5 or float(row[4]) < 0:
                continue  # skip malformed lines
            try:
                value = float(row[4])
            except ValueError:
                continue  # skip if last column isn't a float

            if value == 0.0:
                if current_type is None:
                    # Haven't started a block yet -> hold zeros
                    pending_zeros.append(row)
                else:
                    # Already in a block -> keep zeros inside
                    current_rows.append(row)

            elif value == 1800.0:
                if current_type != "print":
                    # flush previous block
                    if current_type == "reheat":
                        dwell_row = copy.deepcopy(row)
                        dwell_row[4] = 0.0
                        dwell_row[0] = float(row[0]) - epsilon
                        current_rows.append(dwell_row)
                        flush_rows(current_rows, "reheat", reheat_idx)
                        reheat_idx += 1
                    # start new print block with pending zeros first
                    current_rows = pending_zeros[:] + [row]
                    pending_zeros.clear()
                    current_type = "print"
                else:
                    current_rows.append(row)

            else:
                # any other nonzero value -> reheat block
                if current_type != "reheat":
                    # flush previous block
                    if current_type == "print":
                        dwell_row = copy.deepcopy(row)
                        dwell_row[4] = 0.0
                        dwell_row[0] = float(row[0]) - epsilon
                        current_rows.append(dwell_row)
                        flush_rows(current_rows, "print", print_idx)
                        print_idx += 1
                    # start new reheat block with pending zeros first
                    current_rows = pending_zeros[:] + [row]
                    pending_zeros.clear()
                    current_type = "reheat"
                else:
                    current_rows.append(row)

        # Flush last accumulated block
        if current_type == "print":
            flush_rows(current_rows, "print", print_idx)
        elif current_type == "reheat":
            flush_rows(current_rows, "reheat", reheat_idx)


def pick_chunk(values: List[float], layer_idx: int, num_layers: int) -> float:
    if isinstance(values, (int, float)):
        return float(values)  # single value for all layers
    n_chunks = len(values)
    width = num_layers / n_chunks
    cidx = min(int(layer_idx // width), n_chunks - 1)
    return values[cidx]
  
def flatten_and_subtract_off_nominal_dwells(output_dir, d0, d1, num_layers, layer_height, rigid_shift_y_reheat):
    output_dir = Path(output_dir)
    print_dir = output_dir / "print_layers"
    reheat_dir = output_dir / "reheat_layers"

    q = num_layers // 4  # integer floor, stable
    quarter_layers_1based = {q*1, q*2, q*3}

    i = 1
    while (i <= num_layers):
        height = layer_height * float(i)

        nominal_d0 = pick_chunk(d0, i, num_layers-1) 
        nominal_d1 = pick_chunk(d1, i, num_layers-1) 

        if ((i+1) in quarter_layers_1based):
            nominal_d1 = 1200.0

        print_file = print_dir / Path("layer_" + str(i) + "_scan_path.txt")
        with open(print_file, newline="") as f:
            rows = []
            reader = csv.reader(f)
            for row in reader:
                row[3] = height
                rows.append(row)
            
            actual_dwell = float(rows[-1][0]) - float(rows[-2][0])
            print(i, "print: actual", actual_dwell, "nominal", nominal_d0, rows[-1][4], rows[-2][4])
            undwelled_layer_end = float(rows[-2][0]) + (actual_dwell - nominal_d0)
            rows[-1][0] = undwelled_layer_end

        with open(print_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(rows)

        reheat_file = reheat_dir / Path("layer_" + str(i) + "_scan_path.txt")
        with open(reheat_file, newline="") as f:
            rows = []
            reader = csv.reader(f)
            for row in reader:
                row[3] = height
                row[2] = float(row[2]) + rigid_shift_y_reheat
                rows.append(row)
            
            actual_dwell = float(rows[-1][0]) - float(rows[-2][0])
            print(i, "reheat: actual", actual_dwell, "nominal", nominal_d1)
            undwelled_layer_end = float(rows[-2][0]) + (actual_dwell - nominal_d1)
            rows[-1][0] = undwelled_layer_end

        with open(reheat_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerows(rows)

        i = i + 1

if __name__ == "__main__":
    num_layers = 85
    dwell_0 = [111.95093, 46.80592, 98.10942, 48.97092]
    dwell_1 = [55.20886]
    reheat_power = [990.72487, 687.35519, 531.07650, 1340.78436]
    layer_height = 0.00071
    rigid_shift_y_reheat = 0.0158

    print_path = 'print_layers'
    reheat_path = 'reheat_layers'
    
    split_csv("tf_log.csv", "")

    flatten_and_subtract_off_nominal_dwells("", dwell_0, dwell_1, num_layers, layer_height, rigid_shift_y_reheat)

    toolpath_info = get_toolpath_info(print_path, reheat_path, dwell_0, dwell_1, reheat_power,layer_end_time_discretization=5.0)

    toolpath_info['set_dwell_every_n_layers'] = True
    toolpath_info['includes_end_message'] = True

    tpp_clean, layer_end_times, actual_layer_variables = create_toolpath(toolpath_info)

    write_event_series(tpp_clean, toolpath_info['scan_path_out'], toolpath_info['includes_end_message'])
    generate_print_plan_file(toolpath_info, "print_plan.json", actual_layer_variables=actual_layer_variables)
