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

def plot_comparison(prefix, column_index):
    base_path = "./C12"
    best_file = os.path.join(base_path, "best", f"{prefix}.tot")
    not_best_file = os.path.join(base_path, "not_best", f"{prefix}.tot")
    
    if not os.path.exists(best_file) or not os.path.exists(not_best_file):
        print("One or both files do not exist.")
        return
    
    best_data = read_file(best_file)
    not_best_data = read_file(not_best_file)
    
    if best_data.size == 0 or not_best_data.size == 0:
        print("One or both files contain no valid data.")
        return
    
    plt.figure(figsize=(8, 6))
    
    plt.plot(best_data[:, 0], best_data[:, column_index], label=f"Best {column_index}", linestyle='-', color='blue')
    plt.plot(not_best_data[:, 0], not_best_data[:, column_index], label=f"Not Best {column_index}", linestyle='-', color='red')
    
    plt.xlabel("E-out [MeV]")
    plt.ylabel("Cross Section [mb/MeV]")
    plt.title(f"Comparison of {prefix} - Column {column_index}")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 mycode.py <prefix> <column_index>")
    else:
        prefix = sys.argv[1]
        try:
            column_index = int(sys.argv[2])
            plot_comparison(prefix, column_index)
        except ValueError:
            print("Column index must be an integer.")