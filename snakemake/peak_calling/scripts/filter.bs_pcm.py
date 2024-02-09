import argparse
import sys

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Filter lines based on the fourth column.")
    parser.add_argument("-bed", type=str, default="-", help="Path to BED file or '-' for stdin (default: stdin).")
    parser.add_argument("-cell_type", type=str, choices=['mesophyll', 'bundle_sheath'], 
                        help="Specify the cell type for additional filtering.")

    # Parse arguments
    args = parser.parse_args()

    # Process the input
    try:
        if args.bed == "-":
            # If '-' is provided, read from standard input
            process_stream(sys.stdin, args.cell_type)
        else:
            # If a file path is provided, open and process it
            with open(args.bed, 'r') as file:
                process_stream(file, args.cell_type)
    except IOError as e:
        print(f"Error: {e}")

def process_stream(stream, cell_type):
    # Define the allowed strings for the fourth column
    allowed_strings = ["bundle_sheath", "bundle_sheath,unknown",
                       "bundle_sheath,companion_cells_sieve_elements",
                       "bundle_sheath,procambial_meristem", "bundle_sheath,procambium"]

    # Process each line in the input stream
    for line in stream:
        columns = line.strip().split()
        if len(columns) >= 4:
            column_name_split = columns[3].split(';')
            # Apply filtering logic based on cell_type
            if cell_type == 'mesophyll' and column_name_split[1] == 'mesophyll':
                print(line.strip())
            elif cell_type == 'bundle_sheath' and column_name_split[1] in allowed_strings:
                print(line.strip())

if __name__ == "__main__":
    main()






































