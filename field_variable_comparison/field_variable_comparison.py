import numpy as np
import pyvista
import os
import glob

def field_variable_comparison(file1, file2, field_name):
    # We need to loop over each of the MPI domains
    filename_pattern = file1 + '.*.vtu'
    list_of_mpi_rank_files = glob.glob(filename_pattern)
    print("list of files: ", list_of_mpi_rank_files)
    raw_ranks = [int(os.path.basename(f).split('.')[-2]) for f in list_of_mpi_rank_files]
    ranks = sorted(set(raw_ranks))
    print("sorted ranks: ", ranks)
    mpi_domains = max(raw_ranks) + 1

    total_normalized_diff = 0.0

    for i in range(0,mpi_domains):
        full_filename_file1 = file1 + "." + str(i) + ".vtu"
        full_filename_file2 = file2 + "." + str(i) + ".vtu"

        dataset1 = pyvista.read(full_filename_file1)
        dataset1.set_active_scalars("temperature")
        
        dataset2 = pyvista.read(full_filename_file2)
        dataset2.set_active_scalars("temperature")

        # Compute volumes and areas
        sized = dataset1.compute_cell_sizes()

        # Grab volumes for all cells in the mesh
        cell_volumes = sized.cell_data["Volume"]

        celldata1 = dataset1.cells
        celldata2 = dataset2.cells


        # TODO: I need to do this properly -- the easiest way might be interpolating onto the smaller mesh

        # Find the minimum length of the two arrays
        if (len(celldata1) != len(celldata2)):
            print("WARNING: MESHES HAVE DIFFERENT NUMBERS OF CELLS! TRIMMING TO FIT, WHICH MAY INTRODUCE ARTIFACTS")

            min_length = min(len(celldata1), len(celldata2))

            # Trim both arrays to the minimum length
            celldata1 = celldata1[:min_length]
            celldata2 = celldata2[:min_length]

        normalized_diff = np.abs(celldata1 - celldata2)

        # TODO: Need to normalize by cell volume but the values are on the nodes
        #normalized_diff = normalized_diff * cell_volumes[0:len(normalized_diff)-1]

        total_normalized_diff = total_normalized_diff + normalized_diff.sum()

    return total_normalized_diff






