import csv
import sys

def reorder_columns(input_file, output_file, new_order):
    with open(input_file, 'r') as infile:
        reader = list(csv.reader(infile))  # Read all lines into a list

    # Extract headers
    header1 = reader[0]
    header2 = reader[1]
    data_rows = reader[2:]  # Data starts from the third row

    # Remove '#' from the first header row but keep it in the output
    header1_clean = [col.lstrip('#') for col in header1]  # Remove leading '#'

    # Create a mapping of column names to their indexes
    col_index_map = {col: i for i, col in enumerate(header1_clean)}

    # Get the new column order indexes based on user input
    try:
        new_indexes = [col_index_map[col] for col in new_order]
    except KeyError as e:
        print(f"Error: Column '{e.args[0]}' not found in input file.")
        sys.exit(1)

    # Reorder headers and data rows
    reordered_header1 = ['#' + header1_clean[i] for i in new_indexes]  # Re-add '#'
    reordered_header2 = [header2[i] for i in new_indexes]
    reordered_data = [[row[i] for i in new_indexes] for row in data_rows]

    # Write to the output CSV
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(reordered_header1)
        writer.writerow(reordered_header2)
        writer.writerows(reordered_data)

    print(f"Reordered CSV saved as: {output_file}")

# Run the script from command line
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 sortCols.py input.csv output.csv \"COL1,COL2,COL3,...\"")
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    new_order = sys.argv[3].split(",")

    reorder_columns(input_filename, output_filename, new_order)
