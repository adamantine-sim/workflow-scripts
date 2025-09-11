import copy
import numpy as np
import re
import os
import json
import random
from typing import List, Tuple

def get_time_position_power_inp(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    time_position_power = []

    for line in lines:
        if line != '\n':
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

def write_event_series(time_position_power, filename, include_end_message):
    scan_path_string = ""
    
    for entry in time_position_power:
        time = entry[0]
        if time > 1e-4:
            time_string = str(time)
        else:
            time_string = '{:.0e}'.format(time)

        line_string = time_string + ',' + str(entry[1][0]) + ',' + str(entry[1][1]) + ',' + str(entry[1][2]) + ',' + str(entry[2]) + '\n'

        scan_path_string = scan_path_string + line_string

    if include_end_message:
        scan_path_string = scan_path_string + "SCAN_PATH_END"
    else:
        scan_path_string = scan_path_string + "\n"
    
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

def pick_chunk(values: List[float], layer_idx: int, num_layers: int) -> float:
    if isinstance(values, (int, float)):
        return float(values)  # single value for all layers
    n_chunks = len(values)
    width = num_layers / n_chunks
    cidx = min(int(layer_idx // width), n_chunks - 1)
    return values[cidx]
  
def get_toolpath_info(print_path, reheat_path, dwell_0, dwell_1, reheat_power, layer_end_time_discretization=None, layer_limit=None):
    toolpath_info = {}
    toolpath_info['print_path'] = print_path
    toolpath_info['reheat_path'] = reheat_path
    toolpath_info['dwell_0'] = dwell_0 # s
    toolpath_info['dwell_1'] = dwell_1 # s
    toolpath_info['reheat_power'] = reheat_power # W
    toolpath_info['scan_path_out'] = "scan_path.inp"
    toolpath_info['lump_size'] = 2
    toolpath_info['includes_end_message'] = True
    toolpath_info['layer_end_time_discretization'] = layer_end_time_discretization

    toolpath_info['admissible_controls'] = ('reheat_power, dwell_0', 'dwell_1')

    filename_pattern = r'layer_(\d+)_scan_path\.txt'
    filenames_print = get_sorted_layer_files(print_path, filename_pattern)
    filenames_reheat = get_sorted_layer_files(reheat_path, filename_pattern)

    # Apply layer limit if given
    if layer_limit is not None:
        filenames_print = filenames_print[0:layer_limit+1]
        filenames_reheat = filenames_reheat[0:layer_limit+1]

    toolpath_info['num_layers'] = len(filenames_print)

    base_split_layers_print = []
    for file in filenames_print:
        base_split_layers_print.append(get_time_position_power_inp(print_path + '/' + file))
    
    base_split_layers_reheat = []
    for file in filenames_reheat:
        base_split_layers_reheat.append(get_time_position_power_inp(reheat_path + '/' + file))

    # Currently we assume that all layers are the same height
    toolpath_info['layer_height'] = (base_split_layers_print[toolpath_info['lump_size']][1][1][2] - base_split_layers_print[0][1][1][2])/toolpath_info['lump_size']

    toolpath_info['base_split_layers_print'] = base_split_layers_print
    toolpath_info['base_split_layers_reheat'] = base_split_layers_reheat
    
    # By default, select all layers
    toolpath_info['selected_layers'] = (0, toolpath_info['num_layers'])
    return toolpath_info

def generate_print_plan_file(toolpath_info, plan_filename):
    dwell_0       = toolpath_info['dwell_0']
    reheat_power  = toolpath_info['reheat_power']
    dwell_1       = toolpath_info['dwell_1']
    num_layers = toolpath_info.get('num_layers')

    data = []
    for layer_idx in range(num_layers):
        d0 = pick_chunk(dwell_0,      layer_idx, num_layers)
        rp = pick_chunk(reheat_power, layer_idx, num_layers)
        d1 = pick_chunk(dwell_1,      layer_idx, num_layers)

        entry = {
            "layer": layer_idx,
            "power": rp,
            "dwell_0": d0,
            "dwell_1": d1
        }
        data.append(entry)

    with open(plan_filename, "w") as f:
        json.dump(data, f, indent=2)


def create_toolpath(toolpath_info):
    """
    Builds a new scan_path.inp by:
      - loading each peeled layer from print_layers/ and reheat_layers/
      - inserting dwell_0, reheat_power, dwell_1 in equal-width chunks
      - shifting times accordingly
    """
    # 1) Unpack inputs
    print_path                    = toolpath_info['print_path']
    reheat_path                   = toolpath_info['reheat_path']
    dwell_0                       = toolpath_info['dwell_0']
    reheat_power                  = toolpath_info['reheat_power']
    dwell_1                       = toolpath_info['dwell_1']
    includes_end_message          = toolpath_info.get('includes_end_message', True)
    layer_end_time_discretization = toolpath_info.get('layer_end_time_discretization', None)
    set_dwell_every_n_layers      = toolpath_info.get('set_dwell_every_n_layers')
    base_split_layers_print       = toolpath_info.get('base_split_layers_print')
    base_split_layers_reheat      = toolpath_info.get('base_split_layers_reheat')
    num_layers                    = toolpath_info.get('num_layers')
    selected_layers               = toolpath_info.get('selected_layers')  # tuple (start, end)

    # 2) Discover peeled layer files
    filename_pattern = r'layer_(\d+)_scan_path\.txt'
    # grab the raw filenames
    raw_print_files  = get_sorted_layer_files(print_path, filename_pattern)
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
    
    # 3) load the base_split_layers_* if not defined in parameters
    if base_split_layers_print is None:
        base_split_layers_print = [
            get_time_position_power_inp(os.path.join(print_path, fn))
            for fn in filenames_print
        ]
        toolpath_info['base_split_layers_print'] = base_split_layers_print 

    if base_split_layers_reheat is None:
        base_split_layers_reheat = [
            get_time_position_power_inp(os.path.join(reheat_path, fn))
            for fn in filenames_reheat
        ]
        toolpath_info['base_split_layers_reheat'] = base_split_layers_reheat
    
    # 4) Infer the number of layers and selected layers if not set in parameters
    if num_layers is None:
        num_layers = len(base_split_layers_print)
        toolpath_info['num_layers'] = num_layers

    if selected_layers is None:
        selected_layers = (0, num_layers)
        toolpath_info['selected_layers'] = selected_layers
   
    
    # 5) Build the new time–position–power series
    new_tpp = []
    section_start_time = 1e-10
    layer_end_times = []
    # quarters: 1/4, 2/4, 3/4 of the layer count
    q = num_layers // 4  # integer floor, stable
    quarter_layers_1based = {q*1, q*2, q*3}
    
    for layer_idx in range(*toolpath_info['selected_layers']):
        # 5a) Pick parameters for this layer
        d0 = pick_chunk(dwell_0,      layer_idx, num_layers)
        rp = pick_chunk(reheat_power, layer_idx, num_layers)
        d1 = dwell_1 #pick_chunk(dwell_1,      layer_idx, num_layers)
        # If we happen to be at the start of one of the sections, add a 20 min dwell instead of the optimized dwell.
        if set_dwell_every_n_layers and ((layer_idx + 1) in quarter_layers_1based):
            d0 = 1200.0
            
            
        # 5b) Print slice
        layer = base_split_layers_print[layer_idx]
        if layer_idx > 1:
            shifted = shift_time(layer, section_start_time+d0)
        else:
            shifted = shift_time(layer, section_start_time)
        
        # start-of-slice power toggle
        slice_entries = add_power_off_entry(shifted)
        
        # ensure the very last print-point has power=0.0
        if slice_entries:
            t_last, pos_last, _ = slice_entries[-1]
            slice_entries[-1] = (t_last, pos_last, 0.0)    
        new_tpp += slice_entries
        section_start_time = new_tpp[-1][0]
        # 5c) First dwell + Reheat pass + second dwell 
        reheat = base_split_layers_reheat[layer_idx]
        reheat = update_power(reheat, rp)
        reheat = shift_time(reheat, section_start_time)
        new_tpp += reheat

        # End of reheat
        layer_end_time = new_tpp[-1][0]
        pos = new_tpp[-1][1]

        # Dwell 1
        if d1 > 0.0:
            new_tpp += time_position_power_dwell(layer_end_time, pos, d1)
            layer_end_time = new_tpp[-1][0]

        # Snap to discretization (optional)
        if layer_end_time_discretization is not None:
            remainder = layer_end_time % layer_end_time_discretization
            # add only the needed remainder; 0 if already on-grid
            additional = (layer_end_time_discretization - remainder) % layer_end_time_discretization
            if additional > 0.0:
                new_tpp += time_position_power_dwell(layer_end_time, pos, additional)
                layer_end_time = new_tpp[-1][0]

        section_start_time = layer_end_time
        layer_end_times.append(layer_end_time)

    # 6) Clean up and return
    tpp_clean = strip_duplicate_locations(new_tpp)
    #tpp_clean = new_tpp

    return tpp_clean, layer_end_times

def write_toolpath(toolpath_info):
    
    tpp_clean, layer_end_times = create_toolpath(toolpath_info)
    # print("layer end times", layer_end_times)

    write_event_series(tpp_clean, toolpath_info['scan_path_out'], toolpath_info['includes_end_message'])


def generate_control_options(layer_time_discretization, num_forward_sims, sampling_strategy, toolpath_info, planned_controls):
    
    def get_uniform_random_discrete_value(bounds, discrete_step):
        x = random.uniform(bounds[0], bounds[1])
        out = round(x / discrete_step) * discrete_step
        return out


    dwell_0_bounds = (layer_time_discretization, 10*layer_time_discretization)
    dwell_1_bounds = (layer_time_discretization, 10*layer_time_discretization)
    reheat_power_bounds = (0.0, 500.0)

    # The first entry should be the planned control
    modified_toolpath_info_list = []
    modified_toolpath_info_list.append(copy.deepcopy(toolpath_info))
    modified_toolpath_info_list[-1]['dwell_0'] = [planned_controls['dwell_0']]
    modified_toolpath_info_list[-1]['dwell_1'] = [planned_controls['dwell_1']]
    modified_toolpath_info_list[-1]['reheat_power'] = [planned_controls['power']]

    num_variables = 3
    num_adjacent_options = 2 * num_variables + 1


    match sampling_strategy:
        case 'uniform_random':
            for i in range(num_forward_sims-1):   
                
                unique_parameter_set_found = False
                counter = 0

                while not unique_parameter_set_found:
                    # Draw new values
                    dwell_0_val = get_uniform_random_discrete_value(dwell_0_bounds, layer_time_discretization)
                    dwell_1_val = get_uniform_random_discrete_value(dwell_1_bounds, layer_time_discretization)
                    reheat_val = get_uniform_random_discrete_value(reheat_power_bounds, 1.0)

                    # Compare to existing values
                    for entry in modified_toolpath_info_list:
                        if np.isclose(dwell_0_val, entry['dwell_0']) and np.isclose(dwell_1_val, entry['dwell_1']) and np.isclose(reheat_val, entry['reheat_power']):
                            counter = counter + 1
                            if counter > 1000:
                                print('Error: Unable to find enough unique random samples for the control simulations')
                                sys.exit()
                            continue

                    print("Control options:", dwell_0_val, dwell_1_val, reheat_val)

                    modified_toolpath_info_list.append(copy.deepcopy(toolpath_info))
                    modified_toolpath_info_list[-1]['dwell_0'] = [dwell_0_val]
                    modified_toolpath_info_list[-1]['dwell_1'] = [dwell_1_val]
                    modified_toolpath_info_list[-1]['reheat_power'] = [reheat_val]
                    unique_parameter_set_found = True

            return modified_toolpath_info_list
        
        case 'normal_random':
            # TODO
            return modified_toolpath_info_list

        case 'individual_perturbations':
            neighbors_to_sample_from = np.ceil(num_forward_sims/num_adjacent_options)
            # TODO
            return modified_toolpath_info_list
        
        case _:
            print('Error: Invalid sampling strategy chosen:', sampling_strategy)
            sys.exit()
    
    return modified_toolpath_info_list
    
    
if __name__ == "__main__":
    cwd = os.getcwd()
    toolpath_info = {
        'print_path'    : os.path.join(cwd, "print_layers"),
        'reheat_path'   : os.path.join(cwd, "reheat_layers"),
        'dwell_0'       : [20, 30, 40, 50],
        'reheat_power'  : [360, 450, 600, 800],
        'dwell_1'       : 20,
        'scan_path_out' : "scan_path_test.inp",
        'includes_end_message': True,
        'set_dwell_every_n_layers': True
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
