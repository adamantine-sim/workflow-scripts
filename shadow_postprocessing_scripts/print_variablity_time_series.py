import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.interpolate import PchipInterpolator
import csv
import argparse
import os

# NOTE: This script assumes that the "single_point_ensemble_time_series.py" 
# script has been run first to generate the single-print mean temperature CSV files.

# -----------------------------------------------------------------------------
# Setting up the argument parser
# -----------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="Print variability plotter",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--data-directory", default='.', help="path to the directory with the CSV files to be loaded")
parser.add_argument("--filename-pattern", default='mean_p', help="naming pattern for the CSV files to be loaded, integer for the print index and '.csv' extension is assumed.")
args = parser.parse_args()

data_directory = args.data_directory
filename_pattern = args.filename_pattern

# -----------------------------------------------------------------------------
# Plot formatting
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Main script
# -----------------------------------------------------------------------------

plot_filename = 'variability_single_point.png'

# Load the data from the individual prints, not assuming that the time stamps are the same
merged_data = {}

files = os.listdir(data_directory)

max_print_index = -1

for f in files:
    if f.startswith(filename_pattern):
        print_index = int((f.strip(filename_pattern)).strip('.csv'))
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
plt.savefig(plot_filename, format='png')





        
