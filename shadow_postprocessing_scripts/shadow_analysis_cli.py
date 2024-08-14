import argparse
from shadow_analysis import shadow_analysis

parser = argparse.ArgumentParser(description="Digital shadow plotter",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--data-directory", help="location of the files for VisIt to load")
parser.add_argument("-n", "--filename", default='solution', help="base simulation filename")
parser.add_argument("-o", "--output-directory", help="location to write the images to")
parser.add_argument("--previous-print-data-path", help="path to data from previous prints")
parser.add_argument("--sim-field-plot", action='store_true', default=False, help="whether to plot the simulated field")
parser.add_argument("--expt-field-plot", action='store_true', default=False, help="whether to plot the experimental field on the simulation mesh")
parser.add_argument("--single-time-series-plot", action='store_true', default=False, help="whether to plot the time series for a single print")
parser.add_argument("--variability-time-series-plot", action='store_true', default=False, help="whether to plot the time series for a multiple prints")
parser.add_argument("-p", "--point", help="point of interest for time series plots")
parser.add_argument("--visit-path", help="path to the VisIt executable")
parser.add_argument("--print-index", default=0, help='index of the print for variability plots')
parser.add_argument("--rayfile-path", help='the full path to the rayfile used to set the view direction')
args = parser.parse_args()

plot_sim_field = args.sim_field_plot
plot_expt_field = args.expt_field_plot
plot_single_time_series = args.single_time_series_plot
plot_variability_time_series = args.variability_time_series_plot

path_to_adamantine_files = args.data_directory
adamantine_filename = args.filename
output_directory = args.output_directory
previous_print_data_path = args.previous_print_data_path
point_of_interest = args.point
path_to_visit = args.visit_path
print_index = args.print_index
rayfile = args.rayfile_path

print(args)

shadow_analysis(plot_sim_field, plot_expt_field, plot_single_time_series, plot_variability_time_series, path_to_adamantine_files, adamantine_filename, output_directory, previous_print_data_path, point_of_interest, path_to_visit, rayfile, print_index)