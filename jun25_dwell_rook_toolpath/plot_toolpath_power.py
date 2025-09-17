import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
filename_expt = "tf_log.csv"  # <-- replace with your file path
df = pd.read_csv(filename_expt, header=None)  # header=None if no column names in the file

# Select first and last columns
x_expt = df.iloc[:, 0]
y_expt = df.iloc[:, -1]

# Load the CSV file
filename = "scan_path.inp"  # <-- replace with your file path
df = pd.read_csv(filename, header=None)  # header=None if no column names in the file

# Select first and last columns
x = df.iloc[:, 0]
y = df.iloc[:, -1]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(x_expt, y_expt, color='b', marker=".", linestyle="-")
plt.plot(x, y, color='r', marker="o", markerfacecolor='none', linestyle="-")
plt.xlabel("Time")
plt.ylabel("Power")
plt.grid(True)
plt.show()