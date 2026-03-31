import os
import re
import subprocess

def preprocess_stl(inputs):
    ''' creates a copy of the stl file provided from peregrine '''
    ''' rescales and cleans the copied stl file to the template directory '''
    mesh = inputs["mesh"]

    stl_path = mesh["stl_path"]

    stl_file_name = os.path.basename(stl_path)

    # copy the stl to constant/triSurface for openfoam functions
    working_stl_dir = os.path.join(inputs["template"]["template_dir"], "constant", "triSurface")
    working_stl_path = os.path.join(working_stl_dir, stl_file_name)

    os.system(f"mkdir -p {working_stl_dir}")
    os.system(f"cp -rf {stl_path} {working_stl_path}")

    # scale the stl file to meters
    scaling = " ".join(str(mesh["convertToMeters"]) for _ in range(3))
    os.system(f'surfaceTransformPoints "scale=({scaling})" {working_stl_path} {working_stl_path}')

    # generic surface clean (removes ambiguous patches in stl file)
    os.system(f"surfaceClean {working_stl_path} {working_stl_path} 0 0")

    return working_stl_path


def create_background_mesh(inputs, working_stl_path):
    ''' create a background mesh using blockMesh around the stl file '''
    mesh = inputs["mesh"]
    spacing = mesh["spacing"]
    tolerance = mesh["tolerance"]

    # get the bounding box of the stl to create background mesh
    s = subprocess.check_output(f"surfaceCheck {working_stl_path} | grep -i 'Bounding Box :'", shell=True).decode('utf-8')
    bb_str = re.findall('\(([^)]+)', s)

    bb_min = [float(x) - tolerance for x in bb_str[0].split(" ")]
    bb_max = [float(x) + tolerance for x in bb_str[1].split(" ")]
    bb = bb_min + bb_max
    bbDict = {"bb_min": bb_min, "bb_max": bb_max, "bb": bb}

    # set the number of cells in each direction based on the desired spacing
    span = [b - a for (a, b) in zip(bb_min, bb_max)]
    nCells = [round(a / b) for (a, b) in zip(span, spacing)]
    origin = [a + b/2.0 for (a, b) in zip(bb_min, span)]

    # update the background mesh file
    blockMeshDict = os.path.join(inputs["template"]["template_dir"], "system/blockMeshDict")
    lines = open(blockMeshDict, 'r').readlines()

    keys = ["xmin", "ymin", "zmin", "xmax", "ymax", "zmax"]

    for k, key in enumerate(keys):
        for i, line in enumerate(lines):
            if line.startswith(key):
                token = line.replace(";", "").split()[-1]
                lines[i] = line.replace(token, str(bb[k]))

    keys = ["nx", "ny", "nz"]

    for k, key in enumerate(keys):
        for i, line in enumerate(lines):
            if line.startswith(key):
                token = line.replace(";","").split()[-1]                
                lines[i] = line.replace(token, str(nCells[k]))

    with open(blockMeshDict, 'w') as f:
        for line in lines:
            f.write(str(line))

    os.system(f'blockMesh -case {inputs["template"]["template_dir"]}')
    #os.system(f'cd {inputs["template"]["template_dir"]}; blockMesh; cd ..')

    return [origin, bbDict]


def extract_stl_features(inputs, origin):
    ''' extract features from the stl file and set parameters in files '''

    template_dir = inputs["template"]["template_dir"]

    refinement = inputs["mesh"]["refinement"]
    stl_file_name = os.path.basename(inputs["mesh"]["stl_path"])

    # extract the surface features from the stl file
    surfaceFeaturesDict = f"{template_dir}/system/surfaceFeaturesDict"

    os.system(f"foamDictionary -entry surfaces -set '(\"{stl_file_name}\")' {surfaceFeaturesDict}")

    os.system(f"surfaceFeatures -case {template_dir}")

    eMesh_name = stl_file_name.split(".")[0] + ".eMesh"

    # update entries is snappyHexMeshDict
    snappyHexMeshDict = f"{template_dir}/system/snappyHexMeshDict"
    origin = " ".join(list(map(str, origin)))
    print(origin)

    os.system(
        f"foamDictionary -entry geometry/part/file"
        f''' -set '"{stl_file_name}"' {snappyHexMeshDict}''')

    os.system(
        f"foamDictionary -entry castellatedMeshControls/features"
        f''' -set '( {"{"} file "{eMesh_name}"; level {refinement}; {"}"} )' '''
        f"{snappyHexMeshDict}")

    os.system(
        f"foamDictionary -entry castellatedMeshControls/locationInMesh"
        f" -set '( {origin} )' {snappyHexMeshDict}")

    os.system(
        f"foamDictionary -entry castellatedMeshControls/refinementSurfaces/part/level"
        f" -set '( {refinement} {refinement} );' {snappyHexMeshDict}")

    os.system(
        f"foamDictionary -entry castellatedMeshControls/refinementRegions/part/levels"
        f" -set '( ({refinement} {refinement}) );' {snappyHexMeshDict}")

def create_part_mesh(inputs, bbDict):
    ''' create the part mesh '''
    template_dir = inputs["template"]["template_dir"]

    exe = inputs["exe"]

    if exe["nProcs"] <= 1:
        os.system(f"snappyHexMesh -case {template_dir} -overwrite")
    else:
        os.system(
            f"foamDictionary -entry numberOfSubdomains"
            f" -set {exe['nProcs']} {template_dir}/system/decomposeParDict")

        os.system(f"decomposePar -case {template_dir} -force")

        mpi_args = exe['mpiArgs']

        os.system(f"{mpi_args} snappyHexMesh -case {template_dir} -parallel -overwrite")

        os.system(f"reconstructParMesh -case {template_dir} -withZero -constant")

        os.system(f"rm -rf {template_dir}/processor*")

    # move bottom of mesh to z=0 plane
    translation = " ".join(str(t) for t in [0, 0, -bbDict["bb_min"][2]])
    os.system(f"transformPoints -case {template_dir} \"translate=({translation})\"")

    # save a copy of polyMesh for the part in the working directory
    stl_file_name = os.path.basename(inputs["mesh"]["stl_path"])
    polyMesh_copy = f"{template_dir}/{stl_file_name.split('.')[0]}.polyMesh"
    os.system(f"cp -rf {template_dir}/constant/polyMesh {polyMesh_copy}")


def create_cube_mesh(inputs, rve, rve_pad):
    ''' create a background mesh using blockMesh around the stl file '''
    mesh = inputs["mesh"]
    spacing = mesh["spacing"]
    tolerance = mesh["tolerance"]

    # get the bounding box of the stl to create block mesh
    bb_min = [
        rve[0][0] - rve_pad[0], 
        rve[0][1] - rve_pad[1], 
        rve[0][2] - rve_pad[2]
        ] 
    bb_max = [
        rve[1][0] + rve_pad[0], 
        rve[1][1] + rve_pad[1], 
        rve[1][2]
        ] 
    bb = bb_min + bb_max
    bbDict = {"bb_min": bb_min, "bb_max": bb_max, "bb": bb}


    # set the number of cells in each direction based on the desired spacing
    span = [b - a for (a, b) in zip(bb_min, bb_max)]
    nCells = [round(a / b) for (a, b) in zip(span, spacing)]
    origin = [a + b/2.0 for (a, b) in zip(bb_min, span)]

    # update the background mesh file
    blockMeshDict = os.path.join(inputs["template"]["template_dir"], "system/blockMeshDict")
    lines = open(blockMeshDict, 'r').readlines()

    keys = ["xmin", "ymin", "zmin", "xmax", "ymax", "zmax"]

    for k, key in enumerate(keys):
        for i, line in enumerate(lines):
            if line.startswith(key):
                token = line.replace(";", "").split()[-1]
                lines[i] = line.replace(token, str(bb[k]))

    keys = ["nx", "ny", "nz"]

    for k, key in enumerate(keys):
        for i, line in enumerate(lines):
            if line.startswith(key):
                token = line.replace(";","").split()[-1]                
                lines[i] = line.replace(token, str(nCells[k]))

    with open(blockMeshDict, 'w') as f:
        for line in lines:
            f.write(str(line))

    os.system(f'blockMesh -case {inputs["template"]["template_dir"]}')
    #os.system(f'cd {inputs["template"]["template_dir"]}; blockMesh; cd ..')

    return [origin, bbDict]
