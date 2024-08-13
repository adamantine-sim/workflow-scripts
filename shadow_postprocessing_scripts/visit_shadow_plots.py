import sys
import os
import argparse
import csv
import shutil
import time

# Setting up the argument parser
parser = argparse.ArgumentParser(description="Digital shadow plotter",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--data-directory", help="location of the files for VisIt to load")
parser.add_argument("-n", "--filename", default='solution_m0', help="base simulation filename")
parser.add_argument("-r", "--rayfile", help="full path to the rayfile that determines the camera angle")
parser.add_argument("-o", "--output-directory", default='', help="location to write the images to")
parser.add_argument("-t", "--num-iter", help="simulation iteration number")
parser.add_argument("-e", "--experimental-data", action='store_true', help="plots experimental data on the simulation mesh")
args = parser.parse_args()

# Getting some variables out of the argument parser
num_iter = args.num_iter
data_directory = args.data_directory
output_directory = args.output_directory
filename = data_directory + args.filename + '.' + str(num_iter) + '.pvtu'
rayfile_filename = args.rayfile

# Extract the camera position from the rayfile
with open(rayfile_filename, mode='r') as file:
    reader = csv.reader(file)
    next(reader)  # Skip the header line
    second_line = next(reader)  # Read the second line

view_from_rayfile_string = second_line[:3]
view_from_rayfile = [float(x) for x in view_from_rayfile_string]

view = tuple(view_from_rayfile)
view_up = (0, 0, 1)

# Hard-coding the variable for now
variable = 'temperature'

print("VisIt attepting to load " + filename + "...")

OpenDatabase(filename)

print('done.')

output_filename = None
num_plots = None

if (args.experimental_data):
    output_filename = 'experimental_temperature_'
    num_plots = 3

    AddPlot("Pseudocolor", "temperature", 1, 0)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.pointSize = 0.002
    PseudocolorAtts.pointType = PseudocolorAtts.Sphere 
    PseudocolorAtts.renderSurfaces = 0
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 300
    PseudocolorAtts.max = 1000
    PseudocolorAtts.colorTableName = "plasma"
    SetPlotOptions(PseudocolorAtts)

    AddOperator("Threshold", 0)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 1
    ThresholdAtts.boundsInputType = 0
    ThresholdAtts.listedVarNames = ("temperature")
    ThresholdAtts.zonePortions = (1)
    ThresholdAtts.lowerBounds = (-1e37)
    ThresholdAtts.upperBounds = (1e37)
    ThresholdAtts.defaultVarName = "default"
    ThresholdAtts.defaultVarIsScalar = 0
    ThresholdAtts.boundsRange = ("-1e37:1e37")
    SetOperatorOptions(ThresholdAtts, -1, 0)

    AddOperator("Isovolume", 1)
    IsovolumeAtts = IsovolumeAttributes()
    IsovolumeAtts.lbound = 5
    IsovolumeAtts.ubound = 2000
    IsovolumeAtts.variable = "temperature"
    SetOperatorOptions(IsovolumeAtts, 0, 1)

    AddPlot("Subset", "domains", 1, 0)
    SubsetAtts = SubsetAttributes()
    SubsetAtts.colorType = SubsetAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    SubsetAtts.colorTableName = "Default"
    SubsetAtts.invertColorTable = 0
    SubsetAtts.legendFlag = 1
    SubsetAtts.lineWidth = 0
    SubsetAtts.singleColor = (192, 192, 192, 255)
    SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
    SubsetAtts.SetMultiColor(1, (0, 255, 0, 255))
    SubsetAtts.SetMultiColor(2, (0, 0, 255, 255))
    SubsetAtts.subsetNames = ("0", "1", "2")
    SubsetAtts.opacity = 1
    SubsetAtts.wireframe = 0
    SubsetAtts.drawInternal = 0
    SubsetAtts.smoothingLevel = 0
    SubsetAtts.pointSize = 0.05
    SubsetAtts.pointType = SubsetAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    SubsetAtts.pointSizeVarEnabled = 0
    SubsetAtts.pointSizeVar = "default"
    SubsetAtts.pointSizePixels = 2
    SetPlotOptions(SubsetAtts)

    AddOperator("Threshold", 0)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 0
    ThresholdAtts.boundsInputType = 0
    ThresholdAtts.listedVarNames = ("temperature")
    ThresholdAtts.zonePortions = (1)
    ThresholdAtts.lowerBounds = (5)
    ThresholdAtts.upperBounds = (1e+37)
    ThresholdAtts.defaultVarName = "default"
    ThresholdAtts.defaultVarIsScalar = 0
    ThresholdAtts.boundsRange = ("5:1e+37")
    SetOperatorOptions(ThresholdAtts, -1, 0)

    AddPlot("Mesh", "mesh", 1, 0)
    AddOperator("Threshold", 0)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 0
    ThresholdAtts.boundsInputType = 0
    ThresholdAtts.listedVarNames = ("temperature")
    ThresholdAtts.zonePortions = (1)
    ThresholdAtts.lowerBounds = (5)
    ThresholdAtts.upperBounds = (1e+37)
    ThresholdAtts.defaultVarName = "default"
    ThresholdAtts.defaultVarIsScalar = 0
    ThresholdAtts.boundsRange = ("5:1e+37")
    SetOperatorOptions(ThresholdAtts, -1, 0)

    DrawPlots()

else:
    output_filename = 'simulation_temperature_'
    num_plots = 1

    # Plot the temperature
    AddPlot("Pseudocolor", variable)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = 300
    PseudocolorAtts.max = 1000
    PseudocolorAtts.colorTableName = "plasma"
    SetPlotOptions(PseudocolorAtts)


    AddOperator("Threshold", 1)
    ThresholdAtts = ThresholdAttributes()
    ThresholdAtts.outputMeshType = 0
    ThresholdAtts.boundsInputType = 0
    ThresholdAtts.listedVarNames = ("temperature")
    ThresholdAtts.zonePortions = (1)
    ThresholdAtts.lowerBounds = (5)
    ThresholdAtts.upperBounds = (1e+37)
    ThresholdAtts.defaultVarName = "default"
    ThresholdAtts.defaultVarIsScalar = 0
    ThresholdAtts.boundsRange = ("5:1e+37")
    SetOperatorOptions(ThresholdAtts, -1, 1)
    DrawPlots()

# Set the view
View3DAtts = GetView3D()
View3DAtts.viewNormal = view
View3DAtts.viewUp = view_up
SetView3D(View3DAtts)

# Set annotations
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.triadFlag = 1
AnnotationAtts.axes3D.bboxFlag = 0
AnnotationAtts.axes3D.xAxis.title.visible = 1
AnnotationAtts.axes3D.xAxis.title.userTitle = 0
AnnotationAtts.axes3D.xAxis.title.userUnits = 0
AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
AnnotationAtts.axes3D.xAxis.title.units = ""
AnnotationAtts.axes3D.xAxis.label.visible = 1
AnnotationAtts.axes3D.xAxis.label.scaling = 0
AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.xAxis.grid = 0
AnnotationAtts.axes3D.yAxis.title.visible = 1
AnnotationAtts.axes3D.yAxis.title.font.scale = 1
AnnotationAtts.axes3D.yAxis.title.userTitle = 0
AnnotationAtts.axes3D.yAxis.title.userUnits = 0
AnnotationAtts.axes3D.yAxis.title.units = ""
AnnotationAtts.axes3D.yAxis.label.visible = 1
AnnotationAtts.axes3D.yAxis.label.font.scale = 1
AnnotationAtts.axes3D.yAxis.label.scaling = 0
AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.yAxis.grid = 0
AnnotationAtts.axes3D.zAxis.title.visible = 1
AnnotationAtts.axes3D.zAxis.title.userTitle = 0
AnnotationAtts.axes3D.zAxis.title.userUnits = 0
AnnotationAtts.axes3D.zAxis.title.units = ""
AnnotationAtts.axes3D.zAxis.label.visible = 1
AnnotationAtts.axes3D.zAxis.label.font.scale = 1
AnnotationAtts.axes3D.zAxis.label.scaling = 0
AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
AnnotationAtts.axes3D.zAxis.grid = 0
AnnotationAtts.axes3D.setBBoxLocation = 0
AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
AnnotationAtts.axes3D.triadColor = (0, 0, 0)
AnnotationAtts.axes3D.triadLineWidth = 0
AnnotationAtts.axes3D.triadFont = 0
AnnotationAtts.axes3D.triadBold = 1
AnnotationAtts.axes3D.triadItalic = 1
AnnotationAtts.axes3D.triadSetManually = 0
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
AnnotationAtts.databaseInfoTimeScale = 1
AnnotationAtts.databaseInfoTimeOffset = 0
AnnotationAtts.legendInfoFlag = 1
SetAnnotationAttributes(AnnotationAtts)


for plot in range(0, num_plots):
    plotName = GetPlotList().GetPlots(plot).plotName
    legend = GetAnnotationObject(plotName)
    legend.drawMinMax = 0
    legend.drawTitle = 0
    legend.drawLabels = 0
    
    if plot == 0:
        legend.drawLabels = 1
    else:
        # Move additional legends offscreen
        legend.managePosition = 0
        legend.position = (-1000,-1000)    
    
# Save the plot to file
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = str(output_directory)
SaveWindowAtts.fileName = output_filename + str(num_iter)
SaveWindowAtts.family = 0
SaveWindowAtts.format = SaveWindowAtts.PNG  
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.NONE  
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions 
SaveWindowAtts.pixelData = 1

SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

shutil.copyfile(output_directory + output_filename + str(num_iter) + '.png', output_directory + output_filename + 'latest.png')

sys.exit()
