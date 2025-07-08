import sys
import os
import re

def transform_cut(num, val, op, coord, filter_mode=None):
    """
    filter_mode:
        None = processes all files ending in <num>.C
        'exclude_csi' = processes files that DO NOT contain 'csi'
        'include_csi' = processes only files that contain 'csi'
    """
    setpoint_re = re.compile(r'(cutg->SetPoint\(\s*\d+\s*,\s*)([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)(\s*,\s*)([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)(\s*\);)')
    files = [f for f in os.listdir('.') if f.endswith(f'{num}.C')]

    if filter_mode == 'exclude_csi':
        pattern_csi = re.compile(r'^([pdt]_npt_{}|he3_{}|alphas{}|[pdt]_pt_dE_{})\.C$'.format(num, num, num, num))
        files = [f for f in files if pattern_csi.match(f)]
    elif filter_mode == 'include_csi':
        pattern_si = re.compile(r'^(p|d|t)_pt_{}\b|^(p|d|t)nptcsi_{}\b'.format(num, num))
        files = [f for f in files if pattern_si.match(f)]

    if not files:
        print(f"No valid file ending in '{num}.C' found for processing with filter '{filter_mode}'.")
        return

    for filename in files:
        with open(filename, 'r') as f:
            lines = f.readlines()

        new_lines = []
        for line in lines:
            m = setpoint_re.search(line)
            if m:
                x_val = float(m.group(2))
                y_val = float(m.group(4))

                if coord == 'x':
                    x_val = x_val + val if op == 'add' else x_val * val
                elif coord == 'y':
                    y_val = y_val + val if op == 'add' else y_val * val
                else:
                    print(f"Invalid coordinate: {coord}. Use 'x' or 'y'.")
                    return

                new_line = f"{m.group(1)}{x_val:.6g}{m.group(3)}{y_val:.6g}{m.group(5)}\n"
                new_lines.append(new_line)
            else:
                new_lines.append(line)

        with open(filename, 'w') as f:
            f.writelines(new_lines)
        print(f"File '{filename}' updated.")

if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 7:
        print("Usage: transformCut.py <num> <val> <coord> [-a | -m] [-csi | -si]")
        print()
        print("Arguments:")
        print("  <num>      = number in the file suffix (e.g., 16 will match *16.C)")
        print("  <val>      = numeric value to apply")
        print("  <coord>    = 'x' or 'y' (which coordinate to modify)")
        print("  -a         = addition (default if not specified)")
        print("  -m         = multiplication")
        print("  -csi       = exclude files containing 'csi'")
        print("  -si        = exclude si detectors' cuts, i.e. process only only files containing 'csi'")
        sys.exit(1)

    num = sys.argv[1]
    try:
        val = float(sys.argv[2])
    except ValueError:
        print("Invalid value for <val>. Must be a number.")
        sys.exit(1)

    coord = sys.argv[3].lower()
    if coord not in ['x', 'y']:
        print("Invalid coordinate. Use 'x' or 'y'.")
        sys.exit(1)

    # Defaults
    op = 'add'
    filter_mode = None

    for arg in sys.argv[4:]:
        arg = arg.lower()
        if arg == '-a':
            op = 'add'
        elif arg == '-m':
            op = 'multiply'
        elif arg == '-csi':
            filter_mode = 'exclude_csi'
        elif arg == '-si':
            filter_mode = 'include_csi'
        else:
            print(f"Unknown parameter: {arg}")
            sys.exit(1)

    print("Parameters:")
    print(f"  File suffix      : {num}")
    print(f"  Value            : {val}")
    print(f"  Operation        : {'Addition' if op == 'add' else 'Multiplication'}")
    print(f"  Coordinate       : {coord.upper()}")
    print(f"  CSI Filter       : {filter_mode or 'none'}")
    print()

    transform_cut(num, val, op, coord, filter_mode)
