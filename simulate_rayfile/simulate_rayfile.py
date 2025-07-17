import numpy as np
import pyvista
import random
import os
import glob


def simulate_rayfile(vtk_file, ray_filename, field_name, camera_location, num_rays=10):
    '''
    This function creates a rayfile from a simulation result. 
    '''

    # Get surface points and temperatures
    dataset = pyvista.read(vtk_file)
    dataset = dataset.clean(tolerance=1e-6)

    geometry = dataset.extract_geometry()
    surface = geometry.extract_surface()
    surface_points = surface.points
    surface.set_active_scalars(field_name)
    field_data = surface.point_data[field_name]

    # Filter points to keep only points facing the camera with active material
    indices = []
    for i in range(0,len(surface_points)):
        ip, ic = surface.ray_trace(camera_location, surface_points[i])
        if len(ic) < 2 and field_data[i] > 1e-6:
            indices.append(i)

    # Randomly draw from the list
    sampled_indices = random.sample(indices, num_rays)

    rays = []
    for index in sampled_indices:
        ray = (camera_location, surface_points[index], field_data[index])
        rays.append(ray)

    with open(ray_filename, 'w') as rayfile:
        for ray in rays:
            rayfile.write(', '.join(f"{x:.6f}" for x in ray[0]) + ', ' + ', '.join(f"{x:.6f}" for x in ray[1]) + ', ' + str(ray[2]) + '\n')

        rayfile.close()

    # Plotting for debugging
    #pl = pyvista.Plotter()
    #pl.add_mesh(surface, show_edges=True, color='w', opacity=0.5)
    #pl.add_points(surface_points[sampled_indices],render_points_as_spheres=True, color='b', point_size=10)
    #pl.show()




