import numpy as np
import sys
import copy

from ParseROS import parse_ros_gcode

# NOTE: Layers are zero-indexed

def set_start_time_to_zero(time_position_power):
    time_position_power_new = []

    current_start_time = time_position_power[0][0]
    for entry in time_position_power:
        new_time = np.max([entry[0] - current_start_time, 1e-12])
        time_position_power_new.append((new_time, entry[1], entry[2]))

    return time_position_power_new

def trim_initial_power_off(time_position_power):
    time_position_power_trimmed = copy.deepcopy(time_position_power)
    power_off = True
    i = 0
    while power_off:
        if time_position_power[i][2] < 1e-3:
            time_position_power_trimmed.pop(i)
        else:
            power_off = False

        i = i + 1

    return time_position_power_trimmed

def get_time_position_power(file, average_power, time_increment=0.5):

    chunk_time = True

    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    start_time = float(lines[1].split(',')[1])

    time_position_power = []

    last_time = -100.0

    for line in lines[1:]:
        split_line = line.split(',')
        time = float(split_line[1]) - start_time
        time = np.max([time, 1e-12])
        x = float(split_line[2])
        y = float(split_line[3])
        z = float(split_line[4])

        listed_power = float(split_line[9]) * float(split_line[10])

        if listed_power > 1e-3:
            power = average_power
        else:
            power = 0.0

        if chunk_time:
            if (time - last_time) > time_increment:
                time_position_power.append((time, (x,y,z), power))
                last_time = time
            elif power < 1e-3:
                time_position_power.append((time, (x,y,z), power))
                last_time = time
        else:
            time_position_power.append((time, (x,y,z), power))

    return time_position_power

def get_time_position_power_inp(file, average_power, time_increment=0.5):
    chunk_time = True

    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    start_time = float(lines[0].split(',')[0])

    time_position_power = []

    last_time = -100.0

    for line in lines[1:]:
        split_line = line.split(',')
        time = float(split_line[0]) - start_time
        time = np.max([time, 1e-12])
        x = float(split_line[1])
        y = float(split_line[2])
        z = float(split_line[3])

        listed_power = float(split_line[4])

        if listed_power > 1e-3:
            power = average_power
        else:
            power = 0.0

        if chunk_time:
            if (time - last_time) > time_increment:
                time_position_power.append((time, (x,y,z), power))
                last_time = time
            elif power < 1e-3:
                time_position_power.append((time, (x,y,z), power))
                last_time = time
        else:
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

def set_total_layer_time(time_position_power, total_layer_time):
    total_layer_time_set = False

    while not total_layer_time_set:
        entry_time = time_position_power[-1][0]

        if total_layer_time > entry_time:
            time_position_power.append( (total_layer_time, (0.0, 0.0, 0.0), 0.0) )
            total_layer_time_set = True
        else:
            if time_position_power[-1][2] > 1e-3:
                print("ERROR: Requested layer time is too short")
                sys.exit()

            time_position_power.pop()

    return time_position_power

def get_layer_height(time_position_power):
    layer_heights = []
    current_layer_index = 0
    
    for entry in time_position_power:
        if entry[2] > 1e-3:
            z = entry[1][2]

            if len(layer_heights) == 0:
                layer_heights.append(z)

            elif (z - layer_heights[-1] > 1e-3):
                layer_heights.append(z)
                
        if len(layer_heights) == 2:
            continue

    return layer_heights[1] - layer_heights[0]

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
    
    f = open(filename, 'w')
    f.write(scan_path_string)
    f.close()

def set_z_and_shift_time(time_position_power, z, time_shift):
    new_time_position_power = []

    for entry in time_position_power:
        shifted_time = np.max([entry[0] + time_shift, 1e-12])
        shifted_entry = (shifted_time, (entry[1][0], entry[1][1], z), entry[2])
        new_time_position_power.append(shifted_entry)

    return new_time_position_power

def create_composite_scan_path(even_layer, odd_layer, total_height, substrate_height, layer_height, lumped_layers):

    total_layers = int(np.round(total_height/layer_height))

    print("total layer count", total_layers)

    time_position_power = []
    for i in range(0, total_layers):
        if len(time_position_power) > 0:
                last_time = time_position_power[-1][0]
        else:
            last_time = 0.0

        lumped_layer_index = np.ceil( (i+1) / lumped_layers) * lumped_layers

        z_height = substrate_height + lumped_layer_index * layer_height

        first_line = 1
        if i == 0:
            first_line = 0

        if np.mod(i, 2) == 0:
            new_layer = set_z_and_shift_time(even_layer, z_height, last_time)
            time_position_power = time_position_power + new_layer[first_line:]
        else:
             new_layer = set_z_and_shift_time(odd_layer, z_height, last_time)
             time_position_power = time_position_power + new_layer[first_line:]

    return time_position_power

def create_composite_scan_path_arb(layer_list, total_height, substrate_height, layer_height, lumped_layers):

    periodicity = len(layer_list)

    total_layers = int(np.round(total_height/layer_height))

    print("total layer count", total_layers)

    time_position_power = []
    for i in range(0, total_layers):
        if len(time_position_power) > 0:
                last_time = time_position_power[-1][0]
        else:
            last_time = 0.0

        lumped_layer_index = np.ceil( (i+1) / lumped_layers) * lumped_layers

        z_height = substrate_height + lumped_layer_index * layer_height

        period_location = np.mod(i,periodicity)
        new_layer = set_z_and_shift_time(layer_list[period_location], z_height, last_time)

        first_line = 1
        if i == 0:
            first_line = 0
        
        time_position_power = time_position_power + new_layer[first_line:]

    return time_position_power

def create_composite_scan_path_arb_2section(layer_list1, layer_list2, section_switch, total_height, substrate_height, layer_height, lumped_layers):


    total_layers = int(np.round(total_height/layer_height))

    print("total layer count", total_layers)

    time_position_power = []
    for i in range(0, total_layers):
        if len(time_position_power) > 0:
                last_time = time_position_power[-1][0]
        else:
            last_time = 0.0

        if i < section_switch:
            layer_list = layer_list1
        else:
            layer_list = layer_list2

        periodicity = len(layer_list)

        lumped_layer_index = np.ceil( (i+1) / lumped_layers) * lumped_layers

        z_height = substrate_height + lumped_layer_index * layer_height

        period_location = np.mod(i,periodicity)
        new_layer = set_z_and_shift_time(layer_list[period_location], z_height, last_time)

        first_line = 1
        if i == 0:
            first_line = 0
        
        time_position_power = time_position_power + new_layer[first_line:]

    return time_position_power


def main():
    time_increment = 0.25
    wall_height = 0.075

    # Extract stringer layers
    stringer_ros_file = "D2_ROS.csv"
    stringer_average_power = 4500
    stringer_ros_time_position_power = get_time_position_power(stringer_ros_file, stringer_average_power, time_increment)
    layer_height = get_layer_height(stringer_ros_time_position_power)

    print(layer_height)

    stringer_even_time_position_power = extract_single_layer(stringer_ros_time_position_power, 0)
    stringer_odd_time_position_power = extract_single_layer(stringer_ros_time_position_power, 1)

    stringer_even_time_position_power = set_total_layer_time(stringer_even_time_position_power, 60.)
    stringer_odd_time_position_power = set_total_layer_time(stringer_odd_time_position_power, 60.)

    write_event_series(stringer_even_time_position_power, "stringer_even.inp")
    write_event_series(stringer_odd_time_position_power, "stringer_odd.inp")

    # Long stringers (same gcode as regular stringers)
    stringer_long_even_time_position_power = set_total_layer_time(stringer_even_time_position_power, 90.)
    stringer_long_odd_time_position_power = set_total_layer_time(stringer_odd_time_position_power, 90.)

    write_event_series(stringer_long_even_time_position_power, "stringer_long_even.inp")
    write_event_series(stringer_long_odd_time_position_power, "stringer_long_odd.inp")

    # Slow stringers
    stringer_slow_ros_file = "D4_ROS.csv"
    stringer_average_power = 4500
    stringer_slow_ros_time_position_power = get_time_position_power(stringer_slow_ros_file, stringer_average_power, time_increment)
    
    stringer_slow_even_time_position_power = extract_single_layer(stringer_slow_ros_time_position_power, 0)
    stringer_slow_odd_time_position_power = extract_single_layer(stringer_slow_ros_time_position_power, 1)

    stringer_slow_even_time_position_power = set_total_layer_time(stringer_slow_even_time_position_power, 60.)
    stringer_slow_odd_time_position_power = set_total_layer_time(stringer_slow_odd_time_position_power, 60.)

    write_event_series(stringer_slow_even_time_position_power, "stringer_slow_even.inp")
    write_event_series(stringer_slow_odd_time_position_power, "stringer_slow_odd.inp")

    # Extract 2-segment stringer layers
    filename = "D5_gcode.apt"
    rapidSpeed = 5000
    two_segment_average_power = 4500
    dwell = 0.0
    z_offset = 0.012
    discretize_z = True
    z_step = 1.33e-3, 
    two_segment_velocity = 8.0
    event_series_filename = "D5_from_gcode.inp"

    parse_ros_gcode(filename, rapidSpeed, two_segment_average_power, dwell, z_offset, discretize_z, z_step, two_segment_velocity, event_series_filename)

    two_segment_ros_time_position_power = get_time_position_power_inp(event_series_filename, two_segment_average_power, time_increment)

    two_segment_ros_time_position_power = trim_initial_power_off(two_segment_ros_time_position_power)

    two_segment_ros_time_position_power = set_start_time_to_zero(two_segment_ros_time_position_power)

    #layer_height = get_layer_height(two_segment_ros_time_position_power)
    two_segment_even_time_position_power = extract_single_layer(two_segment_ros_time_position_power, 0)
    two_segment_odd_time_position_power = extract_single_layer(two_segment_ros_time_position_power, 1)

    two_segment_even_time_position_power = set_total_layer_time(two_segment_even_time_position_power, 60.)
    two_segment_odd_time_position_power = set_total_layer_time(two_segment_odd_time_position_power, 60.)

    write_event_series(two_segment_even_time_position_power, "two_segment_even.inp")
    write_event_series(two_segment_odd_time_position_power, "two_segment_odd.inp")

    # Extract high curvature stringer layers
    filename = "D7_gcode.apt"
    rapidSpeed = 5000
    high_curvature_average_power = 4500
    dwell = 0.0
    z_offset = 0.012
    discretize_z = True
    z_step = 1.33e-3, 
    high_curvature_velocity = 8.0
    event_series_filename = "D7_from_gcode.inp"

    parse_ros_gcode(filename, rapidSpeed, high_curvature_average_power, dwell, z_offset, discretize_z, z_step, high_curvature_velocity, event_series_filename)

    high_curvature_ros_time_position_power = get_time_position_power_inp(event_series_filename, high_curvature_average_power, time_increment)

    high_curvature_ros_time_position_power = trim_initial_power_off(high_curvature_ros_time_position_power)
    high_curvature_ros_time_position_power = set_start_time_to_zero(high_curvature_ros_time_position_power)

    #layer_height = get_layer_height(two_segment_ros_time_position_power)
    high_curvature_even_time_position_power = extract_single_layer(high_curvature_ros_time_position_power, 0)
    high_curvature_odd_time_position_power = extract_single_layer(high_curvature_ros_time_position_power, 1)

    high_curvature_even_time_position_power = set_total_layer_time(high_curvature_even_time_position_power, 60.)
    high_curvature_odd_time_position_power = set_total_layer_time(high_curvature_odd_time_position_power, 60.)

    write_event_series(high_curvature_even_time_position_power, "high_curvature_even.inp")
    write_event_series(high_curvature_odd_time_position_power, "high_curvature_odd.inp")

    # Weave
    weave_ros_file = "D8_ROS.csv"
    weave_average_power = 4500
    weave_ros_time_position_power = get_time_position_power(weave_ros_file, weave_average_power, time_increment)
    
    weave_even_time_position_power = extract_single_layer(weave_ros_time_position_power, 0)
    weave_odd_time_position_power = extract_single_layer(weave_ros_time_position_power, 1)

    weave_even_time_position_power = set_total_layer_time(weave_even_time_position_power, 60.)
    weave_odd_time_position_power = set_total_layer_time(weave_odd_time_position_power, 60.)

    write_event_series(weave_even_time_position_power, "weave_even.inp")
    write_event_series(weave_odd_time_position_power, "weave_odd.inp")

    # =============================================================================
    # Create the full scan paths
    # =============================================================================
    substrate_height = stringer_ros_time_position_power[0][1][2]
    total_height = wall_height #substrate_height + wall_height

    print("substrate height:", substrate_height)
    print("wall height:", wall_height)
    print("total height:", total_height)

    # Create full D2 scan path
    D2_time_position_power = create_composite_scan_path(stringer_even_time_position_power, stringer_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D2_time_position_power, "D2.inp")

    D2_time_position_power = create_composite_scan_path(stringer_even_time_position_power, stringer_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D2_time_position_power, "D2_lumped.inp")

    # Create full D3 scan path
    D3_time_position_power = create_composite_scan_path(stringer_long_even_time_position_power, stringer_long_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D3_time_position_power, "D3.inp")

    D3_time_position_power = create_composite_scan_path(stringer_long_even_time_position_power, stringer_long_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D3_time_position_power, "D3_lumped.inp")

    # Create full D4 scan path
    D4_time_position_power = create_composite_scan_path(stringer_slow_even_time_position_power, stringer_slow_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D4_time_position_power, "D4.inp")

    D4_time_position_power = create_composite_scan_path(stringer_slow_even_time_position_power, stringer_slow_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D4_time_position_power, "D4_lumped.inp")

    # Create full D5 scan path
    D5_time_position_power = create_composite_scan_path(two_segment_even_time_position_power, two_segment_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D5_time_position_power, "D5.inp")

    D5_time_position_power = create_composite_scan_path(two_segment_even_time_position_power, two_segment_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D5_time_position_power, "D5_lumped.inp")

    # Create full D6 scan path
    layer_list = [stringer_even_time_position_power, two_segment_even_time_position_power, stringer_odd_time_position_power, two_segment_odd_time_position_power]

    D6_time_position_power = create_composite_scan_path_arb(layer_list, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D6_time_position_power, "D6.inp")

    D6_time_position_power = create_composite_scan_path_arb(layer_list, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D6_time_position_power, "D6_lumped.inp")

    # Create full D7 scan path
    D7_time_position_power = create_composite_scan_path(high_curvature_even_time_position_power, high_curvature_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D7_time_position_power, "D7.inp")

    D7_time_position_power = create_composite_scan_path(high_curvature_even_time_position_power, high_curvature_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D7_time_position_power, "D7_lumped.inp")

    # Create full D8 scan path
    D8_time_position_power = create_composite_scan_path(weave_even_time_position_power, weave_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D8_time_position_power, "D8.inp")

    D8_time_position_power = create_composite_scan_path(weave_even_time_position_power, weave_odd_time_position_power, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D8_time_position_power, "D8_lumped.inp")

    # Create full D1 scan path
    layer_list1 = [stringer_even_time_position_power, two_segment_even_time_position_power, stringer_odd_time_position_power, two_segment_odd_time_position_power]

    layer_list2 = [stringer_even_time_position_power, stringer_odd_time_position_power]

    section_switch = 28

    D1_time_position_power = create_composite_scan_path_arb_2section(layer_list1, layer_list2, section_switch, wall_height, substrate_height, layer_height, lumped_layers = 1)
    write_event_series(D1_time_position_power, "D1.inp")

    D1_time_position_power = create_composite_scan_path_arb_2section(layer_list1, layer_list2, section_switch, wall_height, substrate_height, layer_height, lumped_layers = 4)
    write_event_series(D1_time_position_power, "D1_lumped.inp")

if __name__ == "__main__":
    main()