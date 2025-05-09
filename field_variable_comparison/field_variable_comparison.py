import numpy as np
import pyvista
import os
import glob
from scipy.spatial import KDTree

def merge_mpi_decomposition(file1):
    """
    At one point I thought that I needed to merge the data from the individual ranks by hand. It turns out
    that I don't need to. I'm keeping this code here for now, but it can probably be safely deleted.
    """

    # We need to loop over each of the MPI domains
    filename_pattern = file1 + '.*.vtu'
    list_of_mpi_rank_files = glob.glob(filename_pattern)
    #print("list of files: ", list_of_mpi_rank_files)
    raw_ranks = [int(os.path.basename(f).split('.')[-2]) for f in list_of_mpi_rank_files]
    ranks = sorted(set(raw_ranks))
    #print("sorted ranks: ", ranks)
    mpi_domains = max(raw_ranks) + 1

    total_normalized_diff = 0.0

    mpi_subdomains = []

    for i in range(0,mpi_domains):
        full_filename_file1 = file1 + "." + str(i) + ".vtu"
        dataset1 = pyvista.read(full_filename_file1)
        dataset1.set_active_scalars("temperature")
        mpi_subdomains.append(dataset1)

    blocks = pyvista.MultiBlock(mpi_subdomains)
    merged = blocks.combine(merge_points=False)

    #pl = pyvista.Plotter()
    #pl.add_mesh(merged, show_edges=True, cmap='plasma')
    #pl.show()

    return merged

def merged_field_variable_comparison(dataset1, dataset2, field_name):
    # NOTE: There is some difference in the field values for different processor counts. I've confirmed by brute-force search that the issue isn't multiple cells with the same center. Either the simulations depend slightly on the processor count or maybe PyVista is converting to single precision at some point.

    dataset1_cell_centers = dataset1.cell_centers().points
    dataset2_cell_centers = dataset2.cell_centers().points

    tree = KDTree(dataset2_cell_centers)
    distances, indices = tree.query(dataset1_cell_centers)

    total_diff = 0.0

    for i in range(0, len(indices)):
        cell1 = dataset1.get_cell(i)
        cell2 = dataset2.get_cell(indices[i])

        vals1 = dataset1.point_data[field_name][cell1.point_ids]
        vals2 = dataset2.point_data[field_name][cell2.point_ids]

        diff = np.abs(vals1-vals2)

        total_diff = total_diff + diff.sum()
    
    return total_diff

def pvtu_field_variable_comparison(file1, file2, field_name):

    full_filename1 = file1 + ".pvtu"
    dataset1 = pyvista.read(full_filename1)
    dataset1.set_active_scalars(field_name)

    full_filename2 = file2 + ".pvtu"
    dataset2 = pyvista.read(full_filename2)
    dataset2.set_active_scalars(field_name)

    pl = pyvista.Plotter(off_screen=True)
    pl.add_mesh(dataset1, show_edges=True, cmap='plasma')
    #pl.show()
    #pl.save_graphic("field1_out.pdf")

    diff = merged_field_variable_comparison(dataset1, dataset2, field_name)


    return diff    


def field_variable_comparison(file1, file2, field_name):
    
    diff = pvtu_field_variable_comparison(file1, file2, field_name)

    return diff

