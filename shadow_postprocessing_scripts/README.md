# Shadow Postprocessing Scripts
## shadow_postprocessor.py
This plotting script creates up to four types of plots:
1. An image of the simulated temperature for the most recent time step
2. An image of the experimental data projected onto the mesh for the most recent time step
3. The mean temperature and standard deviation for a specified point in a simulation
4. The mean temperature for the most recent print alongside the mean temperature from previous prints

All of the data and some of the plots are obtained using VisIt's CLI utility.

The plotting script takes several arguments, given in the table below.

Command line arguments:
|Argument | Short Version | Functionality|
|---|---------------------|-------------------------|
|--data-directory | -d | Location of the files for VisIt to load|
|--filename | -n | Base filename of the file for VisIt to load, default: 'solution' | 
|--output-directory | -o | Directory to write the plots to | 
|--sim-field-plot |  | Flag to activate the plotting of the simulated temperature field | 
|--expt-field-plot |  | Flag to activate the plotting of the experimental temperature field on the simulation mesh | 
|--single-time-series-plot |  | Flag to activate the plotting of time series of a point for one print | 
|--variability-time-series-plot |  | Flag to activate the plotting of time series of a point compared to previous prints| 
|--point | -p | Point of interest for time series plots | 
|--visit-path | | Path to the VisIt executable | 
|--rayfile-path | | Full path of the rayfile for VisIt to load to determine the view location | 
|--print-index | | Index for the current print (used for comparing the results of different prints) |


Example usage:
```
python3 shadow_postprocessor.py -d /Path/To/Directory/With/Simulation/Results/ --output-directory /Path/Where/You/Want/The/Image/Files/ -n solution -p '"(0.016, 0.00535, 0.0765)"' --visit-path /Path/To/VisIt/visit --print-index 0 --rayfile-filename rayfile.csv --sim-field-plot --expt-field-plot --single-time-series-plot --variability-time-series-plot
```
