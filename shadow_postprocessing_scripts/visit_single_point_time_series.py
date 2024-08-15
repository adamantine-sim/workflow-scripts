import sys
import csv
import argparse
import os
import glob

# --------------------------------------------------------------------------------
# Define the main function for the script
# --------------------------------------------------------------------------------

def write_time_series(filename, max_output, append_file, output_path, skip, variable, point):
    OpenDatabase(filename + " database", 1)

    AddPlot("Pseudocolor", variable)
    DrawPlots()

    time_list = []
    temperature_list = []

    total_time_slider_states = TimeSliderGetNStates()

    start_state = 0
    end_state = total_time_slider_states

    if end_state > max_output:
        end_state = max_output

    # Check the target file to see if it exists and which data point it contains
    if (os.path.exists(output_path)):
        if append_file:
            with open(output_path) as csv_file:
                reader = csv.reader(csv_file)

                num_lines = 0
                for row in reader:
                    num_lines = num_lines + 1

                start_state = num_lines - 1

    print("START_STATE", start_state, "END_STATE", end_state)

    for n in range(start_state, end_state, skip):
        SetTimeSliderState(n)

        Query("Time")
        current_time = GetQueryOutputValue()
        time_list.append(current_time)

        NodePick(point)
        current_temperature = GetPickOutputObject()[variable]
        temperature_list.append(current_temperature)

        if n > max_output:
            break

    rows = zip(time_list, temperature_list)

    if (start_state  > 0):
        with open(output_path, 'a', newline='') as file:
            writer = csv.writer(file)        
            writer.writerows(rows)    

    else:
        with open(output_path, 'w', newline='') as file:
            writer = csv.writer(file)
        
            writer.writerow(['time (s)', 'temperature (K)'])
            writer.writerows(rows)  

    DeleteAllPlots() 

    return 

# --------------------------------------------------------------------------------
# The main script
# --------------------------------------------------------------------------------

# Setting up the argument parser
parser = argparse.ArgumentParser(description="Digital shadow plotter",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--data-directory", help="location of the files for VisIt to load")
parser.add_argument("-n", "--filename", default='solution', help="base simulation filename")
parser.add_argument("-o", "--output-directory", default='', help="location to write the images to")
parser.add_argument("--output-filename", default='time_series', help="filename to write the images to")
parser.add_argument("-p", "--point", help="(x,y,z) point of interest")
parser.add_argument("--skip", default=1, help="outputs to skip")
parser.add_argument("--max-output", default=1e100, help="maximum output to plot to")
parser.add_argument("-a", '--append-existing', action='store_true', default=False, help="if the target CSV filename exists, append to it")
parser.add_argument("--ensemble", action='store_true', default=False, help='whether the script should look for multiple ensemble members')
args = parser.parse_args()

data_directory = args.data_directory
output_directory = args.output_directory
filename = data_directory + args.filename + '.*.pvtu'
point = args.point
output_filename = args.output_filename
max_output = int(args.max_output)
skip = int(args.skip)
append = args.append_existing
ensemble = args.ensemble

base_output_path = output_directory + '/' + output_filename

point = eval(point)

variable = 'temperature'

# Find the ensemble members
if ensemble:
    filename_pattern = data_directory + args.filename + "_m*.*.pvtu"
    list_of_ensemble_files = glob.glob(filename_pattern)
    # Now extract the bounds of the ensemble IDs
    ensemble_id_min = 1e100
    ensemble_id_max = -1e100
    for f in list_of_ensemble_files:
        prefix = data_directory + args.filename + "_m"
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

    for e in range(ensemble_id_min, ensemble_id_max+1):

        filename = data_directory + args.filename + "_m" + str(e) + '.*.pvtu'

        output_path = base_output_path + "_m" + str(e) + ".csv"

        write_time_series(filename, max_output, append, output_path, skip, variable, point)


else:
    filename = data_directory + args.filename + '.*.pvtu'
    output_path = base_output_path + ".csv"
    write_time_series(filename, max_output, append, output_path, skip, variable, point)

sys.exit()
