import sys
import os
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from analysis_17_4_PH.analysis_17_4_PH import objective_17_4_PH
#import analysis_17_4_PH

def main():
    print("Testing 17-4PH analysis...")

    passed = False

    C_to_K = 273.15


    time_series_C = [[0.0, 1800.0+C_to_K, 100.0+C_to_K, 470.0+C_to_K, 475.0+C_to_K,], 
                     [100.0, 1800.0+C_to_K, 160.0+C_to_K, 470.0+C_to_K, 475.0+C_to_K],
                     [100.0, 1800.0+C_to_K, 140.0+C_to_K, 570.0+C_to_K, 475.0+C_to_K],
                     [1800.0+C_to_K, 140.0+C_to_K, 475.0+C_to_K, 475.0+C_to_K, 485.0+C_to_K],
                     [1800.0+C_to_K, 140.0+C_to_K, 560.0+C_to_K, 140.0+C_to_K, 485.0+C_to_K]]

    time_series = np.array(time_series_C)
    time_series = time_series.transpose()

    print(time_series)
    #time_series = np.array([x + C_to_K for x in time_series_C])
   # time_series.reshape((len(time_series_C),2))


    score = objective_17_4_PH(time_series)

    print(score)

    if np.array_equal(score, [2, 0, 0, 3, 1]):
        passed = True
        print("Test passed!")
    else:
        passed = False
        print("Test failed!")


if __name__ == "__main__":
    main()

