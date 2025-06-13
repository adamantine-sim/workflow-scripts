import copy
import numpy as np
import re
import os


def get_time_position_power_inp(file):
    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    time_position_power = []

    for line in lines:
        if line is not '\n':
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
    #print("STRIP")
    out = []
    for entry in time_position_power:
        if len(out) > 0:
            if entry[1] != out[-1][1]:
                out.append(entry)
        else:
            out.append(entry)

    return out


def get_chunked_value(vals, location, chunk_locations):
    val = vals[0]
    for i in range(0,len(chunk_locations)-1):
        if location > chunk_locations[i] and location < chunk_locations[i+1]:
            val = vals[i]
    return val
  
def get_toolpath_info(print_path, reheat_path):
    toolpath_info = {}
    toolpath_info['print_path'] = print_path
    toolpath_info['reheat_path'] = reheat_path
    toolpath_info['dwell_0'] = [10] # s
    toolpath_info['dwell_1'] = [10] # s
    toolpath_info['reheat_power'] = [500] # W
    toolpath_info['scan_path_out'] = "scan_path.inp"
    toolpath_info['lump_size'] = 2
    toolpath_info['includes_end_message'] = True
    
     # Load the individually sliced layers
   
    filename_pattern = 'layer_(\d+)_scan_path\.txt'
    filenames_print = get_sorted_layer_files(print_path, filename_pattern)
    filenames_reheat = get_sorted_layer_files(reheat_path, filename_pattern)

    base_split_layers_print = []
    for file in filenames_print:
        base_split_layers_print.append(get_time_position_power_inp(print_path + '/' + file))
    
    base_split_layers_reheat = []
    for file in filenames_reheat:
        base_split_layers_reheat.append(get_time_position_power_inp(reheat_path + '/' + file))

    # Currently we assume that all layers are the same height
    toolpath_info['layer_height'] = base_split_layers_print[1][1][2] - base_split_layers_print[0][1][2]

    toolpath_info['base_split_layers_print'] = base_split_layers_print
    toolpath_info['base_split_layers_reheat'] = base_split_layers_reheat

    toolpath_info['num_layers'] = len(base_split_layers_print)

    # By default, select all layers (in the control driver sometimes we just need a subset of the layers)
    toolpath_info['selected_layers'] = (0, toolpath_info['num_layers'])

    return toolpath_info

def create_toolpath(toolpath_info):
    dwell_0 = toolpath_info['dwell_0']
    dwell_1 = toolpath_info['dwell_1']
    reheat_power = toolpath_info['reheat_power']
    num_layers = toolpath_info['num_layers']
    
    layer_index = 0

    section_start_time = 1e-10

    num_dwell_0_chunks = len(dwell_0)
    dwell_0_chunk_locations = []
    for i in range(0, num_dwell_0_chunks-1):
        dwell_0_chunk_locations.append(round(i/num_layers * num_layers))

    num_dwell_1_chunks = len(dwell_1)
    dwell_1_chunk_locations = []
    for i in range(0, num_dwell_1_chunks-1):
        dwell_1_chunk_locations.append(round(i/num_layers * num_layers))

    num_reheat_power_chunks = len(reheat_power)
    reheat_power_chunk_locations = []
    for i in range(0, num_reheat_power_chunks-1):
        reheat_power_chunk_locations.append(round(i/num_layers * num_layers))

    new_time_position_power = []

    for layer_index in range(toolpath_info['selected_layers'][0], toolpath_info['selected_layers'][1]):
        layer = toolpath_info['base_split_layers_print'][layer_index]

        dwell_0_this_chunk = get_chunked_value(dwell_0, layer_index, dwell_0_chunk_locations)
        dwell_1_this_chunk = get_chunked_value(dwell_1, layer_index, dwell_1_chunk_locations)
        reheat_power_this_chunk = get_chunked_value(reheat_power, layer_index, reheat_power_chunk_locations)

        # Each layer goes: print, dwell, reheat, dwell

        #print("LAYER", layer_index)

        flat_layer = flatten_layer_z(layer, 10)

        shifted_layer = shift_time(flat_layer, section_start_time)
        #print("shifted layer", shifted_layer)
        layer_to_add = add_power_off_entry(shifted_layer)

        new_time_position_power = new_time_position_power + layer_to_add

        dwell_start_time = new_time_position_power[-1][0]
        dwell_position = new_time_position_power[-1][1]

        first_dwell = time_position_power_dwell(dwell_start_time, dwell_position, dwell_0_this_chunk)
        new_time_position_power = new_time_position_power + first_dwell

        section_start_time = new_time_position_power[-1][0]
        
        # Timings from real print in reheat layers
        if layer_index < num_layers - 1:
            reheat_pass = toolpath_info['base_split_layers_reheat'][layer_index]
            reheat_pass = update_power(reheat_pass, reheat_power_this_chunk)
            reheat_pass = shift_time(reheat_pass, section_start_time)
            new_time_position_power = new_time_position_power + reheat_pass

            dwell_start_time = new_time_position_power[-1][0]
            dwell_position = new_time_position_power[-1][1]

            second_dwell = time_position_power_dwell(dwell_start_time, dwell_position, dwell_1_this_chunk)
            new_time_position_power = new_time_position_power + second_dwell

        # Increment for the next layer
        section_start_time = new_time_position_power[-1][0]
        layer_index = layer_index + 1

    time_position_power = strip_duplicate_locations(new_time_position_power)

    return time_position_power


def write_toolpath(toolpath_info):
    time_position_power = create_toolpath(toolpath_info)
    write_event_series(time_position_power, toolpath_info['scan_path_out'], toolpath_info['includes_end_message'])


if __name__ == "__main__":
    print_path = "print_layers"
    reheat_path = "reheat_layers"
    toolpath_info = get_toolpath_info(print_path, reheat_path)
    write_toolpath(toolpath_info)
