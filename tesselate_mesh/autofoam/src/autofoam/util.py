import json
import pandas as pd
import os
import subprocess

# read JSON formatted input file
def read_inputs(filename):
  try:
    print(f"reading input from {filename}")
    with open(filename) as f:
      data = f.read()
      inputs = json.loads(data)
  except:
    raise Exception(f"Error trying to parse inputs from {filename}")
  return inputs

def convert_peregrine_scanpath(filename, export_path, power=1):
    ''' converts peregrine scan path units to additivefoam scan path units '''
    df = pd.read_csv(filename, delim_whitespace=True)

    # convert X & Y distances to meters
    df["X(m)"] = df["X(mm)"]*1e-3
    df["Y(m)"] = df["Y(mm)"]*1e-3

    # set Z value to zero
    df["Z(m)"] = df["Z(mm)"]*0

    # format columns
    round_cols = ["X(m)", "Y(m)", "Z(m)"]
    df[round_cols] = df[round_cols].round(6)
    for col in round_cols:
        df[col] = df[col].map(lambda x: f'{str(x).ljust(7+len(str(x).split(".")[0]),"0")}')

    # set the laser power
    df["Power(W)"] = df["Pmod"]*power

    # write the converted path to a new file
    df.to_csv(
        export_path, 
        columns=["Mode", "X(m)", "Y(m)", "Z(m)", "Power(W)", "tParam"], 
        sep="\t", 
        index=False)

'''
def check_openfoam_install_status():

    returns status = 0 if no valid OpenFOAM install was found
    and status = 1 if a valid install was found

    install_loc = subprocess.check_output(["which", "blockMesh"])
    if install_loc == "":
    status = 0
    msg = "OpenFOAM (blockMesh) not found in PATH."
    else:
        version = str(install_loc).split(os.path.sep + "OpenFOAM-")[-1].split(os.path.sep)[0]
        msg = f"OpenFOAM version {version} was found."
    if int(version) == 10:
        status = 1
    else:
        status = 0
        msg += " Version 10 is required for preprocessing, generating scripts instead."
    return [status, msg]
'''
