import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from jun25_dwell_rook_toolpath.toolpath_writer import get_toolpath_info
from jun25_dwell_rook_toolpath.toolpath_writer import create_toolpath
from jun25_dwell_rook_toolpath.toolpath_writer import write_toolpath
from jun25_dwell_rook_toolpath.toolpath_writer import generate_print_plan_file
from jun25_dwell_rook_toolpath.toolpath_writer import write_event_series



def write_full_toolpath():
    print("Test: Full toolpath")

    workflow_scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolpath_generator_path = os.path.join(workflow_scripts_path, 'jun25_dwell_rook_toolpath')
    print_path = os.path.join(toolpath_generator_path, 'print_layers')
    reheat_path = os.path.join(toolpath_generator_path, 'reheat_layers')
    dwell_0 = [10]
    dwell_1 = [10]
    reheat_power = [500]

    toolpath_info = get_toolpath_info(print_path, reheat_path, dwell_0, dwell_1, reheat_power)
    write_toolpath(toolpath_info)

    passed = True

    if passed:
        print("Test passed!")
    else:
        print("Test failed!")

    return passed

def write_two_layer_toolpath():
    print("Test: Two-layer toolpath")

    workflow_scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolpath_generator_path = os.path.join(workflow_scripts_path, 'jun25_dwell_rook_toolpath')
    print_path = os.path.join(toolpath_generator_path, 'print_layers')
    reheat_path = os.path.join(toolpath_generator_path, 'reheat_layers')
    dwell_0 = [10]
    dwell_1 = [10]
    reheat_power = [500]
    
    toolpath_info = get_toolpath_info(print_path, reheat_path, dwell_0, dwell_1, reheat_power)
    layer_thickness = toolpath_info['layer_height']
    toolpath_info['selected_layers'] = [0, 2]
    scan_path_filename = "scan_path_2_layer_test.inp"
    template_scan_path_filename = scan_path_filename
    toolpath_info['scan_path_out'] = template_scan_path_filename
    write_toolpath(toolpath_info)

    passed = True

    if passed:
        print("Test passed!")
    else:
        print("Test failed!")

    return passed

def write_four_layer_discretized_toolpath():
    print("Test: Four-layer discretized toolpath")

    workflow_scripts_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    toolpath_generator_path = os.path.join(workflow_scripts_path, 'jun25_dwell_rook_toolpath')
    print_path = os.path.join(toolpath_generator_path, 'print_layers')
    reheat_path = os.path.join(toolpath_generator_path, 'reheat_layers')
    dwell_0 = [10]
    dwell_1 = [10]
    reheat_power = [500]
    
    toolpath_info = get_toolpath_info(print_path, reheat_path, dwell_0, dwell_1, reheat_power, layer_end_time_discretization=12.0)
    layer_thickness = toolpath_info['layer_height']
    toolpath_info['selected_layers'] = [0, 4]
    scan_path_filename = "scan_path_4_layer_discretized_test.inp"
    template_scan_path_filename = scan_path_filename
    toolpath_info['scan_path_out'] = template_scan_path_filename

    tpp_clean, layer_end_times = create_toolpath(toolpath_info)
    write_event_series(tpp_clean, toolpath_info['scan_path_out'], toolpath_info['includes_end_message'])    

    passed = np.allclose(layer_end_times,[98, 242, 338])

    if passed:
        print("Test passed!")
    else:
        print("Test failed!")

    return passed

def write_print_plan():
    print("Test: Write print plan")

    toolpath_info = {
        'dwell_0'       : [20, 30, 40, 50],
        'reheat_power'  : [360, 450, 600, 800],
        'dwell_1'       : [30],
        "num_layers"    : 85
    }

    generate_print_plan_file(toolpath_info, "test_plan.json")

    passed = True

    if passed:
        print("Test passed!")
    else:
        print("Test failed!")

    return passed

def main():
    print("Testing jun25_dwell_rook_toolpath...")

    write_full_toolpath()

    write_two_layer_toolpath()

    write_four_layer_discretized_toolpath()

    write_print_plan()
    


if __name__ == "__main__":
    main()

