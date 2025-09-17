import csv

def find_time_jumps(input_file, threshold=5.0):
    with open(input_file, newline="") as f:
        reader = csv.reader(f)
        prev_time = None
        line_number = 0
        for row in reader:
            line_number += 1
            if len(row) < 1:
                continue
            try:
                t = float(row[0])
            except ValueError:
                continue  # skip bad lines

            if prev_time is not None:
                delta = t - prev_time
                if delta > threshold:
                    print(f"Jump of {delta:.2f} seconds at line {line_number}: {prev_time:.2f} -> {t:.2f}")
            prev_time = t

if __name__ == "__main__":
    find_time_jumps("input.csv", threshold=5.0)
