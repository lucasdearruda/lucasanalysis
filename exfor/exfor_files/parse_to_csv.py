import csv
import re
import sys

def parse_to_csv(input_file, output_file):
    # Read the input file
    with open(input_file, 'r') as file:
        input_text = file.read()

    # Split the input text into lines
    lines = input_text.splitlines()

    # List to store the parsed data
    parsed_data = []

    # Flag to indicate when to start capturing data
    capture_data = False
    headers = []

    # Iterate through each line
    for line in lines:
        # Check if the line indicates the start of the data section
        if line.startswith('#DATA'):
            capture_data = True
            continue
        elif line.startswith('#ENDDATA'):
            capture_data = False
            break

        # If capturing data, process the line
        if capture_data:
            # Remove leading and trailing whitespace
            line = line.strip()

            # Check if the line contains headers
            if line.startswith('#'):
                # Keep the '#' and split the line into columns
                columns = re.split(r'\s+', line)  # Keep the leading '#'
                headers.append(columns)
            else:
                # Split the line into columns based on whitespace
                columns = re.split(r'\s+', line)
                parsed_data.append(columns)

    # Determine the number of numerical columns
    num_columns = len(parsed_data[0]) if parsed_data else 0

    # Extract headers while keeping the `#` symbol
    header1 = headers[0][:num_columns]
    header2 = headers[1][:num_columns]

    # Write the parsed data to the CSV file
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)

        # Write headers with `#` at the start
        writer.writerow(header1)  # First header row
        writer.writerow(header2)  # Second header row

        # Write the parsed data
        writer.writerows(parsed_data)

    print(f"CSV file generated: {output_file}")

# Run from command line
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_to_csv.py input.txt output.csv")
    else:
        input_filename = sys.argv[1]
        output_filename = sys.argv[2]
        parse_to_csv(input_filename, output_filename)
