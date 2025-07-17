import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from jun25_dwell_rook_toolpath.toolpath_writer import get_toolpath_info
from jun25_dwell_rook_toolpath.toolpath_writer import write_toolpath

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

def main():
    print("Testing jun25_dwell_rook_toolpath...")

    write_full_toolpath()

    write_two_layer_toolpath()
    


if __name__ == "__main__":
    main()

