import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulate_rayfile import simulate_rayfile

def main():
    print("Testing simulate_rayfile...")

    vtk_file = 'solution_dwell.63760.pvtu'
    rayfile = 'test_rayfile.csv'
    num_rays = 100
    field_name = 'temperature'
    camera_location = [1.0, 0.0, 0.5]
    simulate_rayfile.simulate_rayfile(vtk_file, rayfile, field_name, camera_location, num_rays)
    
    # Currently we don't have a good way to check the validity of the rayfile
    print("Test passed!")

if __name__ == "__main__":
    main() 