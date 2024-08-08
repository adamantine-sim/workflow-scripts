# Shadow Postprocessing Scripts
## visit_shadow_plots.py
This plotting script creates image files from VisIt plots of the simulated temperature as well as the experimental data projected onto the mesh. Information is passed into the script using command line arguments and the script is run using VisIt's CLI utility.

Command line arguments:
|Argument | Short Version | Functionality|
|---|---------------------|-------------------------|
|--experimental-data | -e | Plots experimental data instead of simulation data |
|--data-directory | -d | Location of the files for VisIt to load|
|--filename | -n | Base filename of the file for VisIt to load, default: 'solution_m0' | 
|--rayfile | -r | Filename of the rayfile for VisIt to load to determine the view location | 
|--num-iter | -i | Simulation iteration number for the file to be loaded |
|--output-directory | | Location to write the image files to | 


Example for plotting simulation results:
```
/Path/To/Visit/visit -cli --no-win -s visit_shadow_plots.py -d /Path/To/Directory/With/Simulation/Results/ --output-directory /Path/Where/You/Want/The/Image/Files/ -t 2000 -r build2_0_0.csv
```

Example for plotting experimental data on the mesh:
```
/Path/To/Visit/visit -cli --no-win -s visit_shadow_plots.py -d /Path/To/Directory/With/Simulation/Results/ --output-directory /Path/Where/You/Want/The/Image/Files/ -t 102928 -n solution.expt -e -r build2_0_0.csv
```