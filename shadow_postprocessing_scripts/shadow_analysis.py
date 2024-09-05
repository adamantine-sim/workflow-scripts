import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.interpolate import PchipInterpolator
import csv
import os
import glob
import shutil
import shlex
import pyvista
import time

def shadow_analysis(plot_sim_field, plot_expt_field, plot_single_time_series, plot_variability_time_series, path_to_adamantine_files, adamantine_filename, output_directory, previous_print_data_path, point_of_interest, path_to_visit, rayfile, print_index):

    # Some hard-coded variables that we may want to open up to users
    scratch_path = '.'
    csv_filename = 'time_series'
    single_print_plot_filename = output_directory + 'single_point_ensemble.png'
    variability_plot_filename = output_directory + 'variability_single_point.png'
    experiment_filename = adamantine_filename + '.expt'
    saved_previous_temperature_prefix = "mean_p"
    max_output_locations = 100
    variable = 'temperature'

    this_file_path = os.path.dirname(os.path.realpath(__file__))

    # ----------------------------------------------------------
    # Plotting parameters
    # ----------------------------------------------------------

    SMALL_SIZE = 22
    MEDIUM_SIZE = 26
    BIGGER_SIZE = 30

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fWidth = 12
    fHeight = 7.5

    # ----------------------------------------------------------
    # Function to get the simulation iteration count
    # ----------------------------------------------------------

    def get_iteration_count(path_prefix):
        filename_pattern = path_prefix + '*.pvtu'
        files = glob.glob(filename_pattern)

        max_iteration_count = -np.inf
        for f in files:
            temp = f.replace('.pvtu', '')
            post_id_index = temp.rfind('.')
            iteration_string = temp[post_id_index+1:]
            iteration_count = int(iteration_string)
            max_iteration_count = max(max_iteration_count, iteration_count)

        return max_iteration_count

    # ----------------------------------------------------------
    # Check if the correct files exist for plots
    # ----------------------------------------------------------
    sim_pattern = path_to_adamantine_files + adamantine_filename + '_m0.*.pvtu'
    files = glob.glob(sim_pattern)
    if len(files) < 2:
        plot_sim_field = False
        plot_single_time_series = False
        plot_variability_time_series = False
        print("No simulation files, skipping all simulation plots.")

    expt_pattern = path_to_adamantine_files + experiment_filename + '.*.pvtu'
    files = glob.glob(expt_pattern)
    if len(files) == 0:
        plot_expt_field = False
        print("No experimental files, skipping experimental temperature field plot.")

    # ----------------------------------------------------------
    # Generating the various types of plots
    # ----------------------------------------------------------

    # Get the time series data from VisIt, if needed
    if (plot_single_time_series or plot_variability_time_series):
        
        data_df = pd.DataFrame()
        time_df = pd.DataFrame()

        # Find the ensemble members
        
        filename_pattern = path_to_adamantine_files + adamantine_filename + "_m*.*.pvtu"
        list_of_ensemble_files = glob.glob(filename_pattern)
        # Now extract the bounds of the ensemble IDs
        ensemble_id_min = 1e100
        ensemble_id_max = -1e100
        for f in list_of_ensemble_files:
            prefix = path_to_adamantine_files + adamantine_filename + "_m"
            temp = f.replace(prefix,'')
            temp2 = temp.replace(".pvtu",'')
            post_id_index = temp2.rfind('.')
            ensemble_id_string = temp2[0:post_id_index]
            ensemble_id = int(ensemble_id_string)
            ensemble_id_min = min(ensemble_id_min, ensemble_id)
            ensemble_id_max = max(ensemble_id_max, ensemble_id)

        if (ensemble_id_min is not 0):
            print("Error: unexpected ensemble members, no 0 member")
            sys.exit()

        filename_ts_pattern = path_to_adamantine_files + adamantine_filename + '_m0.*.pvtu'

        # Get the list of time steps and sort them
        list_of_member_time_step_files = glob.glob(filename_ts_pattern)
        times = [int(os.path.basename(f).split('.')[1]) for f in list_of_member_time_step_files]
        times.sort()

        for e in range(ensemble_id_min, ensemble_id_max+1):

            # Now get the time series
            temperature_list = []
            time_list = []
            for t in times:
                dataset = pyvista.read(path_to_adamantine_files + adamantine_filename + "_m" + str(e) + '.' + str(t) + '.pvtu')
                dataset.set_active_scalars(variable)

                point = eval(eval(point_of_interest))
                dsp = dataset.find_closest_point(point)
                array = dataset.get_array(variable)
                temperature = array[dsp]
                temperature_list.append(temperature)

                # Read the simulated time (not the iteration number) from the VTK files
                f_path = path_to_adamantine_files + adamantine_filename + "_m" + str(e) + '.' + str(t) + '.0.vtu'
                with open(f_path, 'r') as file:
                    for i, line in enumerate(file):
                        if i == 8:  

                            start_tag_end = line.find('>') + 1  # Position right after the opening tag
                            end_tag_start = line.find('</', start_tag_end)  # Position of the closing tag
                            time_entry = float(line[start_tag_end:end_tag_start].strip())
                            time_list.append(time_entry)
                            break

            # Put the data into a pandas dataframe
            data_df['temperature m' + str(e)] = temperature_list

        time_df['time (s)'] = time_list

        std_dev = data_df.std(axis=1)
        mean = data_df.mean(axis=1)    

        mean_filename = output_directory + saved_previous_temperature_prefix + str(print_index) + '.csv'
        with open(mean_filename, 'w', newline='') as file:
            writer = csv.writer(file)
            
            writer.writerow(['time (s)', 'temperature (K)'])
            rows = zip(time_df['time (s)'], mean)

            writer.writerows(rows)   

    # Plot the mean and standard deviation for the temperature of a single print
    if plot_single_time_series:
        fig, ax = plt.subplots()
        spline_model_mean = PchipInterpolator(time_df['time (s)'], mean)
        spline_model_std_dev = PchipInterpolator(time_df['time (s)'], std_dev)


        t_spline = np.linspace(time_df['time (s)'].min(), time_df['time (s)'].max(), 500)
        m_spline = spline_model_mean(t_spline)
        sd_spline = spline_model_std_dev(t_spline)

        coeff = 2.0
        upper_bound_spline = m_spline + coeff * sd_spline
        lower_bound_spline = m_spline - coeff * sd_spline


        plt.plot(t_spline, m_spline, 'k-', label='mean')
        ax.fill_between(t_spline, upper_bound_spline, lower_bound_spline, alpha=0.5, label='2 std. dev.')
        plt.plot(time_df['time (s)'], mean, 'k.')

        plt.legend(ncol=1)

        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (K)")
        fig.set_figwidth(fWidth)
        fig.set_figheight(fHeight)
        plt.savefig(single_print_plot_filename, format='png')

    # Plot the mean of the temperature for the most recent print along with temperature profiles from previous prints
    if plot_variability_time_series:
        merged_data = {}

        files = os.listdir(previous_print_data_path)

        for f in files:
            if f.startswith(saved_previous_temperature_prefix):
                index = int((f.strip(saved_previous_temperature_prefix)).strip('.csv'))
                df = pd.read_csv(previous_print_data_path + f)

                merged_data[index] = df

        # Plot the data
        fig, ax = plt.subplots()

        for key, value in merged_data.items():
            spline_model_mean = PchipInterpolator(value['time (s)'], value['temperature (K)'])
            t_spline = np.linspace(value['time (s)'].min(), value['time (s)'].max(), 500)
            m_spline = spline_model_mean(t_spline)

            color = 'C' + str(key)

            alpha = 0.2

            plt.plot(value['time (s)'], value['temperature (K)'], '.', color=color, alpha=alpha, label=str(key))
            plt.plot(t_spline, m_spline, '-', alpha=alpha, color=color)

        
        spline_model_mean = PchipInterpolator(time_df['time (s)'], mean)
        t_spline = np.linspace(time_df['time (s)'].min(), time_df['time (s)'].max(), 500)
        m_spline = spline_model_mean(t_spline)

        plt.plot(time_df['time (s)'], mean, '.', color='k', alpha=1.0, label='current print')
        plt.plot(t_spline, m_spline, '-', alpha=1.0, color='k')

        #plt.legend(ncol=1)
        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (K)")
        fig.set_figwidth(fWidth)
        fig.set_figheight(fHeight)
        plt.savefig(variability_plot_filename, format='png')

    # Plot the most recent simulation result
    if plot_sim_field:
        # Extract the camera position from the rayfile
        with open(rayfile, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header line
            second_line = next(reader)  # Read the second line

        view_from_rayfile_string = second_line[:3]
        view_from_rayfile = [float(x) for x in view_from_rayfile_string]

        direction_from_rayfile_string = second_line[3:6]
        direction_from_rayfile = [float(x) for x in direction_from_rayfile_string]


        view = tuple(view_from_rayfile)
        view_up = (0, 0, 1)
        direction = tuple(direction_from_rayfile)

        # NEW WAY WITH PYVISTA
        pyvista.start_xvfb()
        pl = pyvista.Plotter()

        # Auto-detect the number of MPI domains
        filename_pattern = path_to_adamantine_files + adamantine_filename + '_m0.*.*.vtu'

        # Get the list of time steps and sort them
        list_of_mpi_ranks = glob.glob(filename_pattern)
        raw_ranks = [int(os.path.basename(f).split('.')[2]) for f in list_of_mpi_ranks]
        ranks = sorted(set(raw_ranks))
        mpi_domains = max(raw_ranks) + 1

        for i in range(0,mpi_domains):
            sim_filename = adamantine_filename + '_m0'
            iteration_number = get_iteration_count(path_to_adamantine_files + sim_filename)

            file_to_plot = path_to_adamantine_files + sim_filename + '.' + str(iteration_number) + '.' + str(i) + '.vtu'

            dataset = pyvista.read(file_to_plot)
            dataset.set_active_scalars("temperature")

            #bounds = dataset.bounds
            #bounding_box_diagonal = np.sqrt( (bounds[0]-bounds[3])**2 + (bounds[1]-bounds[4])**2 + (bounds[2]-bounds[5])**2 )
            bounding_box_diagonal = 0.4
            focal_point = []
            for i in range(0, len(direction)):
                focal_point.append(view[i] + bounding_box_diagonal * direction[i])
            focal_point = tuple(focal_point)

            threshed = dataset.threshold([5, 10000])

            pl.add_mesh(threshed, show_edges=True, cmap='plasma', clim=[0,1700])

        pl.camera.position = view
        pl.camera.viewup = view_up
        pl.camera.focal_point = focal_point

        pl.save_graphic(output_directory + "simulation_temperature_" + str(iteration_number) + ".pdf")
        shutil.copyfile(output_directory + "simulation_temperature_" + str(iteration_number) + ".pdf", output_directory + "simulation_temperature_latest.pdf")  


    # Plot the most recent experimental data on the simulation mesh
    if plot_expt_field:
        
        # Extract the camera position from the rayfile
        with open(rayfile, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header line
            second_line = next(reader)  # Read the second line

        view_from_rayfile_string = second_line[:3]
        view_from_rayfile = [float(x) for x in view_from_rayfile_string]

        direction_from_rayfile_string = second_line[3:6]
        direction_from_rayfile = [float(x) for x in direction_from_rayfile_string]


        view = tuple(view_from_rayfile)
        view_up = (0, 0, 1)
        direction = tuple(direction_from_rayfile)
        
        pyvista.start_xvfb()
        pl = pyvista.Plotter()

    # Get the list of time steps and sort them
        list_of_mpi_ranks = glob.glob(filename_pattern)
        raw_ranks = [int(os.path.basename(f).split('.')[2]) for f in list_of_mpi_ranks]
        ranks = sorted(set(raw_ranks))
        mpi_domains = max(raw_ranks) + 1

        for i in range(0,mpi_domains):
            sim_filename = adamantine_filename + '.expt'
            iteration_number = get_iteration_count(path_to_adamantine_files + sim_filename)

            file_to_plot = path_to_adamantine_files + sim_filename + '.' + str(iteration_number) + '.' + str(i) + '.vtu'

            dataset = pyvista.read(file_to_plot)

            #bounds = dataset.bounds
            #bounding_box_diagonal = np.sqrt( (bounds[0]-bounds[3])**2 + (bounds[1]-bounds[4])**2 + (bounds[2]-bounds[5])**2 )
            bounding_box_diagonal = 0.4
            focal_point = []
            for i in range(0, len(direction)):
                focal_point.append(view[i] + bounding_box_diagonal * direction[i])
            focal_point = tuple(focal_point)

            # This shows the mesh
            pl.add_mesh(dataset, show_edges=True, color='w')

            dataset.set_active_scalars("temperature")

            # Now I want to find points of interest
            ps = dataset.cast_to_pointset()
            tps = ps.threshold([0, 100000])
            if (len(tps.points) > 0):
                pl.add_points(tps, render_points_as_spheres=True, point_size=10.0, cmap='plasma', clim=[0,1700])

        pl.camera.position = view
        pl.camera.viewup = view_up
        pl.camera.focal_point = focal_point

        pl.save_graphic(output_directory + "experimental_temperature_" + str(iteration_number) + ".pdf")
        shutil.copyfile(output_directory + "experimental_temperature_" + str(iteration_number) + ".pdf", output_directory + "experimental_temperature_latest.pdf")  



