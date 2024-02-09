from sys import argv
import argparse
import sys
import os



def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass






def get_parser():
    parser = argparse.ArgumentParser(description='Given a list of markers \
    prepares  them to be used in the R script annotate_cells.v3.R.')

    parser.add_argument('-o', "--output",  help="base name to be used in the script", 
            required=True, dest='o')
    parser.add_argument('-i', "--input", \
            help="File of markers to split",  \
            required=True, dest='i')
    parser.add_argument('-hd', "--header", \
        action='store_true', dest='hd')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    
    header_line = ["chr","star","end","geneID","name", "type"]
    if args.hd == True:
        cleaned_file_lines = [header_line]
    else:
        cleaned_file_lines = []
    
    with open(args.i, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split("\t")
            take_base_lines = cleaned_line[0:6]
            cleaned_file_lines.append(take_base_lines)

    remove_if_exists((args.o))
    with open(args.o, "a+") as f:
        for line in cleaned_file_lines:
            f.write("\t".join(line))
            f.write("\n")




