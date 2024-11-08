#------------------------------------------------------------------------#
#             Gcode to Event Series Converter for Hybrid AM, Okuma       #
#                    Written by: Rylan Paye, ORNL Intern, 2021           #
#                    Modified by: Yousub Lee 2021/10/18                  #
#                    Modified by: Ashley Gannon 2024/06/05               # 
#                      Oak Ridge National Laboratory                     # 
#------------------------------------------------------------------------#

import re
import numpy as np
import matplotlib.pyplot as plt

# These should be user inputs ------------------------------------------
filename= "ExampleOkuma.MIN" 
beadHeight = 00.76 #1.40    # mm
beadWidth  = 3.5 #6.0     # mm
printSpeed = 900    # mm/min
rapidSpeed = 925.2331 # mm/min
scale = 1e-3; # Adamantine needs conversion to m
# -----------------------------------------------------------------------

gcode=open(filename,'r')
gcodeLines=gcode.readlines()

xdata = np.zeros(len(gcodeLines))
ydata = np.zeros(len(gcodeLines))
zdata = np.zeros(len(gcodeLines))
gdata = np.zeros(len(gcodeLines))
fdata = np.zeros(len(gcodeLines)) # Laser speed control
edata = np.zeros(len(gcodeLines)) # Laser power on/off
pdata = np.zeros(len(gcodeLines)) # Laser power control
mdata = np.zeros(len(gcodeLines)) # Wire on/off
start = np.zeros(len(gcodeLines))
end = np.zeros(len(gcodeLines))

i = 0
for line in gcodeLines:
    xdataStr = re.findall(r'X(\-)?(\d+)?(\.)?(\d+)?', line.strip())
    ydataStr = re.findall(r'Y(\-)?(\d+)?(\.)?(\d+)?', line.strip())
    zdataStr = re.findall(r'Z(\-)?(\d+)?(\.)?(\d+)?', line.strip())  
    fdataStr = re.findall(r'F(\d+)',line.strip())
    lPdataStr = re.search(r'/LPWW=(\d+)',line.strip())
    if lPdataStr:
        laserPower = lPdataStr.group(1)
    pdataStr = re.findall(r'/LPW=(\d+)', line.strip())
    ppdataStr = re.findall(r'/LPW=LPWW', line.strip())
    gdataStr = re.findall(r'G(\d+)',line.strip())   
    startStr = re.findall('5X_begin.txt',line.strip()) 
    endStr   = re.findall(r'files_x\\job_end',line.strip())
    # Set the feedrate
    if fdataStr:
       printSpeed = float(fdataStr[-1])
    # Set x, y, and z if there is no update, save the last value
    if xdataStr:
        if xdataStr[0][1] or xdataStr[0][3]:
            xdataStr = xdataStr[0][0]+xdataStr[0][1]+xdataStr[0][2]+xdataStr[0][3]
            xdata[i] = float(xdataStr)
    else:
        xdata[i] = xdata[i-1]
    if ydataStr:
        if ydataStr[0][1] or ydataStr[0][3]:
            ydataStr = ydataStr[0][0]+ydataStr[0][1]+ydataStr[0][2]+ydataStr[0][3]
            ydata[i] = float(ydataStr)
    else:
        ydata[i] = ydata[i-1]
    if zdataStr:
        if zdataStr[0][1] or zdataStr[0][3]:
            zdataStr = zdataStr[0][0]+zdataStr[0][1]+zdataStr[0][2]+zdataStr[0][3]
            zdata[i] = float(zdataStr)
    else:
        zdata[i] = zdata[i-1]
    # G command
    if gdataStr:
        gdata[i] = gdataStr[0]
    if gdataStr:
        if gdataStr[0]=='1':
            start[i] = i
    if endStr:
        end[i] = i
    # Update counter before power since power ind is offset by 1 for adamantine
    i+=1
    # if /LPW = 0, power is 0
    if len(pdataStr)!=0:
        pdata[i-1] = 0
    # if /LPW = LPWW grab the laserPower
    elif len(ppdataStr)!=0:
        pdata[i-1] = laserPower    
    # if the power doesn't change, keep the old one
    else:
        pdata[i-1] = pdata[i-2]
    # if the laser is on, you are printing
    if pdata[i-1] > 0:
        fdata[i-1] = printSpeed
    # if the laser is off, you are traveling
    else:
        fdata[i-1] = rapidSpeed

# Filter out everything that doesn't happen during the print
start = np.where(start!=0)[0][0]
end   = np.where(end!=0)[0][0]
gcode.close() # Done with this
x = xdata[start:end]
y = ydata[start:end]
z = zdata[start:end]
e = edata[start:end]
p = pdata[start:end]
f = fdata[start:end]
g = gdata[start:end]

# calculate segment times
t = np.zeros(len(x))
mask = np.ones(len(t),dtype=bool)
for k in range(1,len(x)):
    dx = x[k]-x[k-1]
    dy = y[k]-y[k-1]
    dz = z[k]-z[k-1]
    
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    vel = f[k]/60
    dt = dist/vel
    t[k] = t[k-1]+dt
    if k != len(x)-1:
        dz_next = z[k+1]-z[k]
    else:
        dz_next = 0
    if dt == 0: # and dz_next == 0:
        mask[k] = False

x = x[mask]
y = y[mask]
z = z[mask]
e = e[mask]
p = p[mask]
f = f[mask]
t = t[mask]
g = g[mask] 
e = np.roll(e,0)

printTime = np.zeros(len(z))
for ii in range(0,len(e)-2):
    dx = x[ii+1]-x[ii]
    dy = y[ii+1]-y[ii]
    dz = z[ii+1]-z[ii]
    
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    vel = f[ii]/60
    dt = dist/vel
    if p[ii]==10:
        p[ii] = 0
    else:
        printTime[ii+1] = printTime[ii] + dt

p = np.roll(p, -1)  # Shift power values back by one index
p[-1] = 0  # Set the last power value to 0
printTime[0] = 1e-12 # Set first time to small number that is not 0

writeFile = filename.split('.')
writeFile = writeFile[0] + "_EventSeries.inp"
file = open(writeFile,'w')

for j in range(0,len(printTime)-1):
    file.write("{:9.4E},{:.6f},{:.6f},{:.5f},{:.1f}\n".format(printTime[j],x[j]*scale,y[j]*scale,z[j]*scale,p[j]))
    if j==len(printTime)-2:
        print("Data written to file.")

file.close()

