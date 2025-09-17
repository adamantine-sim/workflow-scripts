import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
filename = "print_layers/layer_1_scan_path.txt"  # <-- replace with your file path
df = pd.read_csv(filename, header=None)  # header=None if no column names in the file

x = df.iloc[:, 1]
y = df.iloc[:, 2]

# Load the CSV file
filename = "reheat_layers/layer_1_scan_path.txt"  # <-- replace with your file path
df = pd.read_csv(filename, header=None)  # header=None if no column names in the file

x_rh = df.iloc[:, 1]
y_rh = df.iloc[:, 2]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(x_rh, y_rh, color='b', marker="o", markerfacecolor='none', linestyle="-")
plt.plot(x, y, color='r', marker="o", markerfacecolor='none', linestyle="-")
plt.xlabel("Time")
plt.ylabel("Power")
plt.grid(True)
plt.show()