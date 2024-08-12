import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import subprocess
from scipy.interpolate import PchipInterpolator
import csv

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

plot_filename = 'single_point_ensemble.png'

point_of_interest = '"(0.016, 0.00535, 0.0765)"'

print_index = 0

# VisIt inputs
path_to_visit = '/Applications/VisIt.app/Contents/Resources/bin/visit'
path_to_adamantine_files = '/Users/71d/Documents/workspace/DockerWorkspace_adamantine/aug23_sns_da/'
scratch_path = '.'

adamantine_base_filename = 'solution_m'
ensemble_size = 3


# Extract the data from the VTK files
data = pd.DataFrame()
time = pd.DataFrame()

for member in range(0, ensemble_size):
    adamantine_filename = adamantine_base_filename + str(member)

    csv_filename = 'time_series_m' + str(member) + '.csv'

    visit_command = path_to_visit + ' -cli --no-win -s visit_single_point_time_series.py -d ' + path_to_adamantine_files + ' --output-directory ' + scratch_path + ' --output-filename ' + csv_filename + ' -p ' + point_of_interest + ' -n ' + adamantine_filename + ' --max-output ' + str(10)

    print("VISIT COMMAND", visit_command)

    subprocess.run(visit_command, shell=True)

    df = pd.read_csv(csv_filename)
    field_name = "temperature m" + str(member)
    data[field_name] = df['temperature (K)']
    time['time (s)'] = df['time (s)']


# Calculate summary statistics
std_dev = data.std(axis=1)
mean = data.mean(axis=1)

# Write the mean to file
mean_filename = 'mean_p' + str(print_index) + '.csv'
with open(mean_filename, 'w', newline='') as file:
    writer = csv.writer(file)
    
    writer.writerow(['time (s)', 'temperature (K)'])
    rows = zip(time['time (s)'], mean)

    writer.writerows(rows)    

# Plot the mean and standard deviation for one print
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
plt.savefig(plot_filename, format='png')
