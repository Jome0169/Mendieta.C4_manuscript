from sys import argv
import argparse
import sys
import os

def remove_file(arg1):
    """TODO: Docstring for remove_file(.

    :arg1: TODO
    :returns: TODO

    """
    try:
        os.remove(arg1)
    except OSError:
        pass


def read_file(file_name, n_value):
    """TODO: Docstring for read_file.
    :returns: TODO

    """
    counter = 0
    kept_lines = []
    with open(file_name, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split()
            if cleaned_line[0] == "chr":
                pass
            elif  cleaned_line[0] != "chr":
                if counter < int(n_value):
                    kept_lines.append(cleaned_line)
                    counter += 1
                if counter == n_value:
                    pass
    return(kept_lines)


def grab_header(file_name):
    """Reads in the first line of the first file

    :file_name: TODO
    :returns: TODO

    """
    header_names  = []
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith("chr"):
                cleaned_line = line.strip().split()
                header_names.append(cleaned_line)
            else:
                pass
    return(header_names)

     


def get_parser():
    parser = argparse.ArgumentParser(description='Given a list of markers \
    prepares  them to be used in the R script annotate_cells.v3.R.')

    parser.add_argument('-DA', "--DA_files",  help="Differential Accessability Files", 
            required=True, dest='DA', nargs ="+")
    parser.add_argument('-n', "--number",  help="Number of DA genes to keep per file", 
        required=True, dest='n', type = int)
    parser.add_argument('-o', "--output", help="Output File",  \
            required=False, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    header_values = grab_header(args.DA[0])
    kept_n_lines = [read_file(i, args.n) for i in args.DA]

    def flatten(t):
        return [item for sublist in t for item in sublist]
    flattened = flatten(kept_n_lines)

    if args.o == None:
        print('\t'.join(header_values[0]))
        for gene in flattened:
            print('\t'.join(gene))

    elif args.o != None:
        remove_file(args.o)
        with open(args.o, 'a+') as f:
            f.write("\t".join(header_values[0]))
            f.write('\n')
            for gene in flattened:
                f.write('\t'.join(gene))
                f.write('\n')
