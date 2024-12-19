import re
import numpy as np

#------------------------------------------------------------------------#
#             Gcode to Event Series Converter for Hybrid AM, IMTS        #
#              Written by: Ashley Gannon & Mithulan Paramanathan & Steve DeWitt #
#                      Oak Ridge National Laboratory                     # 
#------------------------------------------------------------------------#

def parse_ros_gcode(filename, rapidSpeed, laserPower, dwell, z_offset, discretize_z, z_step, velocity, event_series_filename):

    feedrate = velocity * 60.

    gcode = open(filename,'r')
    gcodeLines = gcode.readlines()

    x = 12345678*np.ones(len(gcodeLines)) 
    y = 12345678*np.ones(len(gcodeLines))
    z = 12345678*np.ones(len(gcodeLines))
    f = 12345678*np.ones(len(gcodeLines)) # laser speed
    p = 12345678*np.ones(len(gcodeLines)) # laser power
    t = 12345678*np.ones(len(gcodeLines)) # time elapsed for each segment

    i_vals = 12345678*np.ones(len(gcodeLines))
    j_vals = 12345678*np.ones(len(gcodeLines))
    k_vals = 12345678*np.ones(len(gcodeLines))

    last_rapid = False
    deposit = False

    i = 0

    for line in gcodeLines:
        """
        if line.strip().startswith("FEDRAT"):
            match = re.search(r'MMPM\s*,\s*(\d+)', line)
            if match:
                feedrate = float(match.group(1))
        """
            
        if line.strip().startswith("RAPID"):
            # laser power should be 0 for rapid movements. We need to track the rapid on the line prior to GOTO
            last_rapid = True
        if line.strip().startswith("CALSUB/START_DEPO"):
            deposit = True
        if line.strip().startswith("CALSUB/STOP_DEPO"):
            deposit = False
            if i != 1:
                p[i] = 0  
        if line.strip().startswith("GOTO") and deposit == True:
            # Using regular expression to find coordinates after '/'
            match = re.search(r'/\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)', line)
            if match:
                x[i] = float(match.group(1))
                y[i] = float(match.group(2))
                z[i] = float(match.group(3))
                i_vals[i] = float(match.group(4))
                j_vals[i] = float(match.group(5))
                k_vals[i] = float(match.group(6))
            f[i] = feedrate
            if i != 0:
                dx = x[i] - x[i-1]
                dy = y[i] - y[i-1]
                dz = z[i] - z[i-1]
                
                dist = np.sqrt(dx**2 + dy**2 + dz**2)
                vel = f[i]/60
                dt = dist/vel
                t[i] = t[i-1] + dt
            else:
                t[i] = 5e-12 #a small non-zero number
            p[i] = laserPower
            last_rapid = False
            i += 1
        elif line.strip().startswith("GOTO") and last_rapid == True:
            # Using regular expression to find coordinates after '/'
            match = re.search(r'/\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)', line)
            if match:
                x[i] = float(match.group(1))
                y[i] = float(match.group(2))
                z[i] = float(match.group(3))
                i_vals[i] = float(match.group(4))
                j_vals[i] = float(match.group(5))
                k_vals[i] = float(match.group(6))
            f[i] = rapidSpeed
            if i == 0:
                t[i] = 5e-12 #a small non-zero number
            else:
                dx = x[i] - x[i-1]
                dy = y[i] - y[i-1]
                dz = z[i] - z[i-1]
                
                dist = np.sqrt(dx**2 + dy**2 + dz**2)
                vel = f[i]/60
                dt = dist/vel
                t[i] = t[i-1] + dt
            p[i] = 0 # laser is off for rapid moves
            last_rapid = False
            i += 1
        elif line.strip().startswith("GOTO") and deposit == False:
            # Using regular expression to find coordinates after '/'
            match = re.search(r'/\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)\s*,\s*([-+]?\d*\.\d+|\d+)', line)
            if match:
                x[i] = float(match.group(1))
                y[i] = float(match.group(2))
                z[i] = float(match.group(3))
                i_vals[i] = float(match.group(4))
                j_vals[i] = float(match.group(5))
                k_vals[i] = float(match.group(6))
            f[i] = feedrate
            if i == 0:
                t[i] = 5e-12 #a small non-zero number
            else:
                dx = x[i] - x[i-1]
                dy = y[i] - y[i-1]
                dz = z[i] - z[i-1]
                
                dist = np.sqrt(dx**2 + dy**2 + dz**2)
                vel = f[i]/60
                dt = dist/vel
                t[i] = t[i-1] + dt
            p[i] = 0 # laser is off for rapid moves
            last_rapid = False
            i += 1

    x = x[x != 12345678]
    y = y[y != 12345678]
    z = z[z != 12345678]
    t = t[t != 12345678]
    p = p[p != 12345678]

    # Identify the axis of material addition using i,j,k values
    # NOTE: this only works when fully aligned with an axis (i,j,or k = 1)
    power_on_indices = np.where(p != 0)[0]
    i_vals_on = i_vals[power_on_indices]
    j_vals_on = j_vals[power_on_indices]
    k_vals_on = k_vals[power_on_indices]

    # If i or j is ±1, the object needs to be aligned with the z-axis
    if np.any(np.isin(i_vals_on, [-1])):
        # aligned in the -x direction, rotate 90 degrees clockwise about x
        x, z = z, -x
    elif np.any(np.isin(i_vals_on, [1])):
        # aligned in the +x direction, rotate 270 degrees clockwise about x 
        x, z = -z, x
    elif np.any(np.isin(j_vals_on, [1])):
        # aligned in the +y direction, rotate 270 degrees clockwise about y
        y, z = -z, y
    elif np.any(np.isin(j_vals_on, [-1])):
        # aligned in the -y direction, rotate 90 degrees clockwise about y
        y, z = z, -y
        
    '''# Print in the positive z direction
    if np.any(np.diff(z) < 0):
        # If z is decreasing, flip the arrays
        x = x[::-1]
        y = y[::-1]
        z = z[::-1]
        f = f[::-1]
        p = p[::-1]
        p = np.roll(p, -1)  # Shift power values back by one index
        t = t[-1] - t[::-1]  # Adjust the time array'''

    if discretize_z:
        for idx in range(0, len(z)):
            z[idx] = np.ceil(z[idx]/z_step) * z_step

    t[0] = np.max([float(t[0]), 1e-12])
        
    p = np.roll(p,-1)
    p[-1] = 0  # Set the last power value to 0   
    writeFile = event_series_filename
    file = open(writeFile,'w')
    for j in range(0,len(t)):
        file.write("{:.4e},{:.6f},{:.6f},{:.5f},{:.1f}\n".format(t[j],x[j]*1e-3,y[j]*1e-3,z[j]*1e-3+z_offset,p[j]))
        if j==len(t)-2:
            print("Data written to file.")

if __name__ == "__main__":

    filename = "D5_gcode.apt"

    rapidSpeed = 500.0   # mm/min Input from ROS settings is 6000. Doesn't hit that.
    laserPower = 4.456E+3    # W	
    dwell = 0.01           # s
    z_offset = 1.49e-3
    discretize_z = True
    z_step = 1.49e-3 * 2.
    velocity = 8.0 # Feedrate from code overridden by ROS setting (mm/min)

    event_series_filename = "D5_from_gcode.inp"

    parse_ros_gcode(filename, rapidSpeed, laserPower, dwell, z_offset, discretize_z, z_step, velocity, event_series_filename)
