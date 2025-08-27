import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_file(filepath):
    data = []
    with open(filepath, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                try:
                    values = list(map(float, line.split()))
                    data.append(values)
                except ValueError:
                    continue
    return np.array(data)

def plot_comparison(folder, prefix, column_index):
    base_path = folder  # Now the base path is dynamic, using the provided folder
    best_files = [f for f in os.listdir(os.path.join(base_path, "best")) if f.startswith(prefix) and f.endswith(".tot")]
    not_best_files = [f for f in os.listdir(os.path.join(base_path, "not_best")) if f.startswith(prefix) and f.endswith(".tot")]

    if not best_files or not not_best_files:
        print(f"No matching files found for {prefix}")
        return

    # Create the output directory if it doesn't exist
    output_dir = os.path.join("comparison", folder)
    os.makedirs(output_dir, exist_ok=True)

    for best_file, not_best_file in zip(best_files, not_best_files):
        best_file_path = os.path.join(base_path, "best", best_file)
        not_best_file_path = os.path.join(base_path, "not_best", not_best_file)
        
        if not os.path.exists(best_file_path) or not os.path.exists(not_best_file_path):
            print(f"File not found: {best_file_path} or {not_best_file_path}")
            continue
        
        best_data = read_file(best_file_path)
        not_best_data = read_file(not_best_file_path)
        
        if best_data.size == 0 or not_best_data.size == 0:
            print(f"Invalid data in file: {best_file_path} or {not_best_file_path}")
            continue
        
        plt.figure(figsize=(8, 6))
        plt.plot(best_data[:, 0], best_data[:, column_index], label=f"Best {best_file}", linestyle='-', color='blue')
        plt.plot(not_best_data[:, 0], not_best_data[:, column_index], label=f"Not Best {not_best_file}", linestyle='-', color='red')

        plt.xlabel("E-out [MeV]")
        plt.ylabel("Cross Section [mb/MeV]")
        plt.title(f"Comparison of {prefix} - {best_file} vs {not_best_file} - Column {column_index}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        # Save each plot as a PNG file in the comparison/C12/ folder
        output_file = os.path.join(output_dir, f"{prefix}_{os.path.basename(best_file)}_{os.path.basename(not_best_file)}_comparison.png")
        plt.savefig(output_file)
        plt.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 mycode.py <folder> <prefix> <column_index>")
    else:
        folder = sys.argv[1]
        prefix = sys.argv[2]
        try:
            column_index = int(sys.argv[3])
            plot_comparison(folder, prefix, column_index)
        except ValueError:
            print("Column index must be an integer.")
