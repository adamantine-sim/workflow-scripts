import os
import shutil

from melt_pool_analysis.melt_pool_analysis import melt_pool_analysis

def main():
    print("Testing melt-pool-analysis...")
    this_file_path = os.path.dirname(os.path.realpath(__file__))

    path_to_adamantine_files = this_file_path + "/nov24_A2_reference_solution/"
    adamantine_filename = "solution"

    output_directory = this_file_path + "/scratch"

    if os.path.isdir(output_directory) == True:
        shutil.rmtree(output_directory)

    os.mkdir(output_directory)

    melt_pool_analysis(path_to_adamantine_files, adamantine_filename, output_directory)


if __name__ == "__main__":
    main()

