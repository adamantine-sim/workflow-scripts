import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.interpolate import PchipInterpolator
import csv
import os
import glob

# ----------------------------------------------------------
# Handling input arguments
# ---------------------------------------------------------- 
parser = argparse.ArgumentParser(description="Digital shadow plotter",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--data-directory", help="location of the files for VisIt to load")
parser.add_argument("-n", "--filename", default='solution', help="base simulation filename")
parser.add_argument("-o", "--output-directory", help="location to write the images to")
parser.add_argument("--sim-field-plot", action='store_true', default=False, help="whether to plot the simulated field")
parser.add_argument("--expt-field-plot", action='store_true', default=False, help="whether to plot the experimental field on the simulation mesh")
parser.add_argument("--single-time-series-plot", action='store_true', default=False, help="whether to plot the time series for a single print")
parser.add_argument("--variability-time-series-plot", action='store_true', default=False, help="whether to plot the time series for a multiple prints")
parser.add_argument("-p", "--point", help="point of interest for time series plots")
parser.add_argument("--visit-path", help="path to the VisIt executable")
parser.add_argument("--print-index", default=0, help='index of the print for variability plots')
parser.add_argument("--rayfile-filename", help='name of the rayfile, assumed to be in the data-directory')
args = parser.parse_args()

plot_sim_field = args.sim_field_plot
plot_expt_field = args.expt_field_plot
plot_single_time_series = args.single_time_series_plot
plot_variability_time_series = args.variability_time_series_plot

path_to_adamantine_files = args.data_directory
adamantine_filename = args.filename
output_directory = args.output_directory
point_of_interest = args.point
path_to_visit = args.visit_path
print_index = args.print_index
rayfile = args.rayfile_filename


# Some hard-coded variables that we may want to open up to users
scratch_path = '.'
csv_filename = 'time_series'
single_print_plot_filename = output_directory + 'single_point_ensemble.png'
variability_plot_filename = output_directory + 'variability_single_point.png'
saved_temperature_profile_filename = "mean_p"
experiment_filename = adamantine_filename + '.expt'
max_output_locations = 100


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
    filename_pattern = path_prefix + '.*.pvtu'
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
sim_pattern = path_to_adamantine_files + adamantine_filename + '.*.pvtu'
files = glob.glob(sim_pattern)
print(files)
if len(files) == 0:
    plot_sim_field = False
    plot_single_time_series = False
    plot_variability_time_series = False
    print("No simulation files, skipping all simulation plots.")

expt_pattern = path_to_adamantine_files + adamantine_filename + '.expt.*.pvtu'
files = glob.glob(expt_pattern)
print(files)
if len(files) == 0:
    plot_expt_field = False
    print("No experimental files, skipping experimental temperature field plot.")

# ----------------------------------------------------------
# Generating the various types of plots
# ----------------------------------------------------------

# Get the time series data from VisIt, if needed
if (plot_single_time_series or plot_variability_time_series):
    visit_command = path_to_visit + ' -cli -nowin -s visit_single_point_time_series.py -d ' + path_to_adamantine_files + ' --output-directory ' + scratch_path + ' --output-filename ' + csv_filename + ' -p ' + point_of_interest + ' -n ' + adamantine_filename + ' --max-output ' + str(max_output_locations) + ' -a --ensemble'

    print("VISIT COMMAND", visit_command)

    subprocess.run(visit_command, shell=True)

    data = pd.DataFrame()
    time = pd.DataFrame()

    # Determine the number of ensemble members
    n = 0
    num_members = 0 
    while n < 100:
        found_member_csv = os.path.exists(csv_filename + '_m' + str(n) + '.csv')
        if not found_member_csv:
            num_members = n
            break
        n = n + 1

    for member in range(0, num_members):
        df = pd.read_csv(csv_filename + '_m' + str(member) + '.csv')
        field_name = "temperature m" + str(member)
        data[field_name] = df['temperature (K)']
        time['time (s)'] = df['time (s)']

    std_dev = data.std(axis=1)
    mean = data.mean(axis=1)    

    mean_filename = 'mean_p' + str(print_index) + '.csv'
    with open(mean_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        
        writer.writerow(['time (s)', 'temperature (K)'])
        rows = zip(time['time (s)'], mean)

        writer.writerows(rows)    

# Plot the mean and standard deviation for the temperature of a single print
if plot_single_time_series:
    fig, ax = plt.subplots()
    spline_model_mean = PchipInterpolator(time['time (s)'], mean)
    spline_model_std_dev = PchipInterpolator(time['time (s)'], std_dev)


    t_spline = np.linspace(time['time (s)'].min(), time['time (s)'].max(), 500)
    m_spline = spline_model_mean(t_spline)
    sd_spline = spline_model_std_dev(t_spline)

    coeff = 2.0
    upper_bound_spline = m_spline + coeff * sd_spline
    lower_bound_spline = m_spline - coeff * sd_spline


    plt.plot(t_spline, m_spline, 'k-', label='mean')
    ax.fill_between(t_spline, upper_bound_spline, lower_bound_spline, alpha=0.5, label='2 std. dev.')
    plt.plot(time['time (s)'], mean, 'k.')

    plt.legend(ncol=1)

    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    fig.set_figwidth(fWidth)
    fig.set_figheight(fHeight)
    plt.savefig(single_print_plot_filename, format='png')

# Plot the mean of the temperature for the most recent print along with temperature profiles from previous prints
if plot_variability_time_series:
    merged_data = {}

    files = os.listdir(scratch_path)

    max_print_index = -1

    for f in files:
        if f.startswith(saved_temperature_profile_filename):
            print_index = int((f.strip(saved_temperature_profile_filename)).strip('.csv'))
            df = pd.read_csv(f)

            merged_data[print_index] = df

            if max_print_index < print_index:
                max_print_index = print_index


    # Plot the data
    fig, ax = plt.subplots()

    for key, value in merged_data.items():
        spline_model_mean = PchipInterpolator(value['time (s)'], value['temperature (K)'])
        t_spline = np.linspace(value['time (s)'].min(), value['time (s)'].max(), 500)
        m_spline = spline_model_mean(t_spline)

        color = 'C' + str(key)

        alpha = 0.2
        if key == max_print_index:
            alpha = 1.0
            color = 'k'

        plt.plot(value['time (s)'], value['temperature (K)'], '.', color=color, alpha=alpha, label=str(key))
        plt.plot(t_spline, m_spline, '-', alpha=alpha, color=color)


    #plt.legend(ncol=1)

    plt.xlabel("Time (s)")
    plt.ylabel("Temperature (K)")
    fig.set_figwidth(fWidth)
    fig.set_figheight(fHeight)
    plt.savefig(variability_plot_filename, format='png')

# Plot the most recent simulation result
if plot_sim_field:
    sim_filename = adamantine_filename + '_m0'

    iteration_number = get_iteration_count(path_to_adamantine_files + sim_filename)

    visit_command = path_to_visit + ' -cli -nowin -s visit_shadow_plots.py -d ' + path_to_adamantine_files + ' --output-directory ' + output_directory + " -t " + str(iteration_number) + ' -r ' + rayfile + ' -n ' + sim_filename

    print("VISIT COMMAND", visit_command)

    subprocess.run(visit_command, shell=True)

# Plot the most recent experimental data on the simulation mesh
if plot_expt_field:
    iteration_number = get_iteration_count(path_to_adamantine_files + experiment_filename)

    visit_command = path_to_visit + ' -cli -nowin -s visit_shadow_plots.py -d ' + path_to_adamantine_files + ' --output-directory ' + output_directory + " -t " + str(iteration_number) + ' -r ' + rayfile + ' -n ' + experiment_filename + ' -e '

    print("VISIT COMMAND", visit_command)

    subprocess.run(visit_command, shell=True)


