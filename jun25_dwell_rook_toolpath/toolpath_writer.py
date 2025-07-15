import copy
import numpy as np
import re
import os
from typing import List, Tuple

def get_time_position_power_inp(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    time_position_power = []

    for line in lines:
        split_line = line.split(',')
        time = float(split_line[0])
        x = float(split_line[1])
        y = float(split_line[2])
        z = float(split_line[3])
        power = float(split_line[4])
        time_position_power.append((time, (x,y,z), power))

    return time_position_power

def extract_single_layer(time_position_power, target_layer):
    current_layer_height = None
    current_layer_index = 0
    current_layer_start_time = 0.0
    extracted_layer = []
    
    for entry in time_position_power:
        if entry[2] > 1e-3:
            z = entry[1][2]

            if current_layer_height == None:
                current_layer_height = z

            if (z - current_layer_height > 1e-3):
                current_layer_height = z
                current_layer_index = current_layer_index + 1
                current_layer_start_time = entry[0]

        if current_layer_index == target_layer:
            shifted_time = np.max([entry[0] - current_layer_start_time, 1e-12])
            shifted_entry = (shifted_time, (entry[1][0], entry[1][1], 0.0), entry[2])
            extracted_layer.append(shifted_entry)

    return extracted_layer

def write_event_series(time_position_power, filename):
    scan_path_string = ""
    for entry in time_position_power:
        time = entry[0]
        if time > 1e-4:
            time_string = str(time)
        else:
            time_string = '{:.0e}'.format(time)

        line_string = time_string + ',' + str(entry[1][0]) + ',' + str(entry[1][1]) + ',' + str(entry[1][2]) + ',' + str(entry[2]) + '\n'

        scan_path_string = scan_path_string + line_string

    scan_path_string = scan_path_string + "SCAN_PATH_END"
    
    f = open(filename, 'w')
    f.write(scan_path_string)
    f.close()

def split_into_layers(time_position_power, epsilon=0.001):
    # Assumes that the only z changes are between layers, although robot motions could have z moves too
    split_layers = []
    z_layer_start = time_position_power[0][1][2]

    layer = []
    for entry in time_position_power:
        z = entry[1][2]
        if z > z_layer_start+epsilon:
            z_layer_start = z
            split_layers.append(layer)
            layer = [entry]
        else:
            layer.append(entry)

    return split_layers   


def time_position_power_dwell(start_time, position, duration):

    #print("start time", start_time)
    #print("duration", duration)
    dwell_segment = [(start_time+duration, position, 0.0)]

    return dwell_segment

def update_power(time_position_power, power):
    out = []

    for entry in time_position_power:
        #print("entry", entry)
        new_entry = (entry[0], entry[1], power)
        out.append(new_entry)

    return out

def shift_time(layer, new_start_time, eps=1e-12):
    current_start_time = layer[0][0]
    shift = new_start_time - current_start_time + 1e-12

    new_layer = []
    for entry in layer:
        new_layer.append((entry[0]+shift, entry[1], entry[2]))

    return new_layer

def add_power_off_entry(layer, eps=1e-10):
    new_layer = [(layer[0][0], layer[0][1], 0.0)]
    new_layer = new_layer + [(layer[0][0]+ eps, layer[0][1], layer[0][2])]
    new_layer = new_layer + layer[1:]

    return new_layer

def get_sorted_layer_files(directory, file_pattern):
    pattern = re.compile(r'layer_(\d+)_scan_path\.txt')
    matched_files = []

    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            layer_number = int(match.group(1))
            matched_files.append((layer_number, filename))

    # Sort by the extracted number
    matched_files.sort(key=lambda x: x[0])

    # Return only the filenames in order
    return [filename for _, filename in matched_files]


def flatten_layer_z(layer, reference_index=0):
    z_value = layer[reference_index][1][2]

    new_layer = []
    for entry in layer:
        new_layer.append((entry[0], (entry[1][0], entry[1][1], z_value), entry[2]))

    return new_layer

def strip_duplicate_locations(time_position_power):
    out = []
    for entry in reversed(time_position_power):
        if len(out) > 0:
            if entry[1] != out[-1][1]:
                out.append(entry)
        else:
            out.append(entry)

    return list(reversed(out))


def get_chunked_value(vals, location, chunk_locations):
    val = vals[0]
    for i in range(0,len(chunk_locations)-1):
        if location > chunk_locations[i] and location < chunk_locations[i+1]:
            val = vals[i]
    return val
  
def get_toolpath_info(toolpath_path):
    toolpath_info = {}
    toolpath_info['toolpath_path'] = toolpath_path
    toolpath_info['dwell_0'] = [10] # s
    toolpath_info['dwell_1'] = [10] # s
    toolpath_info['reheat_power'] = [500] # W
    toolpath_info['scan_path_out'] = "scan_path.inp"
    toolpath_info['lump_size'] = 2
    # Load the individually sliced layers without dwells or reheat passes
   
    filename_pattern = 'layer_(\d+)_scan_path\.txt'
    filenames = get_sorted_layer_files(toolpath_path, filename_pattern)

    base_split_layers = []
    for file in filenames:
        base_split_layers.append(get_time_position_power_inp(toolpath_path + '/' + file))

    # Currently we assume that all layers are the same height
    toolpath_info['layer_height'] = base_split_layers[1][1][2] - base_split_layers[0][1][2]

    toolpath_info['base_split_layers'] = base_split_layers
    
    toolpath_info['num_layers'] = len(base_split_layers)

    # By default, select all layers
    toolpath_info['selected_layers'] = (0, toolpath_info['num_layers'])
    return toolpath_info



def write_toolpath(toolpath_info):
    """
    Builds a new scan_path.inp by:
      - loading each peeled layer from print_layers/ and reheat_layers/
      - inserting dwell_0, reheat_power, dwell_1 in equal-width chunks
      - shifting times accordingly
      - writing out the final event series
    """
    # 1) Unpack inputs
    print_path    = toolpath_info['print_path']
    reheat_path   = toolpath_info['reheat_path']
    dwell_0       = toolpath_info['dwell_0']
    reheat_power  = toolpath_info['reheat_power']
    dwell_1       = toolpath_info['dwell_1']

    # 2) Discover peeled layer files
    filename_pattern = r'layer_(\d+)_scan_path\.txt'
    #filenames_print  = get_sorted_layer_files(print_path, filename_pattern)
    #filenames_reheat = get_sorted_layer_files(reheat_path, filename_pattern)
    
    # grab the raw filenames
    raw_print_files = get_sorted_layer_files(print_path, filename_pattern)
    raw_reheat_files = get_sorted_layer_files(reheat_path, filename_pattern)
    # then sort them by the integer captured in the filename
    filenames_print = sorted(
        raw_print_files,
        key=lambda fn: int(re.match(r'layer_(\d+)_scan_path\.txt', fn).group(1))
    )
    filenames_reheat = sorted(
        raw_reheat_files,
        key=lambda fn: int(re.match(r'layer_(\d+)_scan_path\.txt', fn).group(1))
    )
    
    # 3) Load the time/position/power for each peeled layer
    base_split_layers_print = [
        get_time_position_power_inp(os.path.join(print_path, fn))
        for fn in filenames_print
    ]
    
    base_split_layers_reheat = [
        get_time_position_power_inp(os.path.join(reheat_path, fn))
        for fn in filenames_reheat
    ]
   
    # 4) Record metadata
    num_layers = len(base_split_layers_print)
    toolpath_info['num_layers'] = num_layers
    toolpath_info['selected_layers'] = (0, num_layers)
    # layer height (assumes uniform layering)
    h0 = base_split_layers_print[0][1][2]
    h1 = base_split_layers_print[1][1][2]
    toolpath_info['layer_height'] = h1 - h0
    toolpath_info['base_split_layers_print'] = base_split_layers_print
    toolpath_info['base_split_layers_reheat'] = base_split_layers_reheat

    # 5) Helper: pick the correct chunk value for a given layer
    def pick_chunk(values: List[float], layer_idx: int) -> float:
        n_chunks = len(values)
        width = num_layers / n_chunks
        cidx = min(int(layer_idx // width), n_chunks - 1)
        return values[cidx]

    # 6) Build the new time–position–power series
    new_tpp = []
    section_start_time = 1e-10
    for layer_idx in range(num_layers):
        # 6a) Pick parameters for this layer
        d0 = pick_chunk(dwell_0,      layer_idx)
        rp = pick_chunk(reheat_power, layer_idx)
        d1 = pick_chunk(dwell_1,      layer_idx)

        # 6b) Print slice
        layer = base_split_layers_print[layer_idx]
        shifted = shift_time(layer, section_start_time)
        # start-of-slice power toggle
        slice_entries = add_power_off_entry(shifted)
        # ensure the very last print-point has power=0.0
        if slice_entries:
            t_last, pos_last, _ = slice_entries[-1]
            slice_entries[-1] = (t_last, pos_last, 0.0)    
        new_tpp += slice_entries
        section_start_time = new_tpp[-1][0] + d0
        
        # 6c) Skip dwell/reheat after the final printed layer NOTE: NO DATA, CHANGE EXP & REMOVE
        if layer_idx == num_layers - 1:
            continue

        # 6d) First dwell + Reheat pass + second dwell (skip after last layer)
        reheat = base_split_layers_reheat[layer_idx]
        reheat = update_power(reheat, rp)
        reheat = shift_time(reheat, section_start_time)
        
        new_tpp += reheat

        t_d1 = new_tpp[-1][0]
        pos  = new_tpp[-1][1]
        new_tpp += time_position_power_dwell(t_d1, pos, d1)
        section_start_time = new_tpp[-1][0]
    # 7) Clean up and write out
    
    tpp_clean = strip_duplicate_locations(new_tpp)
    write_event_series(tpp_clean, toolpath_info['scan_path_out'])

if __name__ == "__main__":
    """
    Quick test harness for write_toolpath():
      - 4 chunks of dwell_0 = [20, 30, 40, 50]
      - 4 chunks of reheat_power = [360, 450, 600, 800]
      - 1 chunk of dwell_1 = [30]
    Outputs scan_path_test.inp in this folder.
    """
    # assume this script lives next to print_layers/ and reheat_layers/
    cwd = os.getcwd()
    toolpath_info = {
        'print_path'    : os.path.join(cwd, "print_layers"),
        'reheat_path'   : os.path.join(cwd, "reheat_layers"),
        'dwell_0'       : [20, 30, 40, 50],
        'reheat_power'  : [360, 450, 600, 800],
        'dwell_1'       : [30],
        'scan_path_out' : "scan_path_test.inp"
    }

    print("Running write_toolpath test with:")
    print(f"  dwell_0       = {toolpath_info['dwell_0']}")
    print(f"  reheat_power  = {toolpath_info['reheat_power']}")
    print(f"  dwell_1       = {toolpath_info['dwell_1']}")
    print(f"  writing -> {toolpath_info['scan_path_out']}")

    # call the core function
    write_toolpath(toolpath_info)

    if os.path.exists(toolpath_info['scan_path_out']):
        print(":) scan_path_test.inp successfully written")
    else:
        print(":( failed to write scan_path_test.inp")
