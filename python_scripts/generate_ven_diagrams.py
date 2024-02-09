import matplotlib
matplotlib.use('Agg')
import pybedtools
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import argparse
import pybedtools
import argparse
import os
import itertools
import copy


def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass



def read_in_bed_file(bed_file):
    """
    Reads in a bed file using the pybedtools standard appraoch.
    """
    def intersecting_feature(feature, index):
        """TODO: Docstring for intersecting_feature(feature, L.
        :returns: TODO
        """
        return(feature[index]) != '.'

    try:
        read_file = pybedtools.BedTool(bed_file)
    except: 
        print("File DOES NOT EXISTS")
        exit(1)

    return read_file 


def shared_peaks(bed_file_1, bed_file_2):
    """
    :returns: TODO

    """
    intersection = bed_file_1 + bed_file_2
    middle_number = intersection.count()
    return middle_number


def unique_peaks(bed_file_1, bed_file_2):
    """TODO: Docstring for unique_peaks.

    :arg1: TODO
    :returns: TODO

    """
    bed_1_only = bed_file_1 - bed_file_2
    bed_2_only = bed_file_2 - bed_file_1

    count_1 = bed_1_only.count()
    count_2 = bed_2_only.count()
    
    return (count_1, count_2)



def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
            replicate peak calls, as well as unqiue peaks to each replicate and \
            outputs said peaks. ')
    parser.add_argument('-bed1','--bed_file1', help='bed1 file',\
            nargs="+", required=True, dest='bed')

    parser.add_argument('-header_name','--headers', help='bed3 file', \
          nargs="+", required=False, dest='heads', type=str)

    parser.add_argument('-title','--title', help='bed3 file', \
            required=False, dest='title', type=str)

    parser.add_argument('-o','--output_name', help='output', \
            required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Read in tissue names
    tissue_list = args.heads

    #Read Files
    bed_list = [read_in_bed_file(i) for i in list(args.bed)]

    bed_1 = bed_list[0]
    bed_2 = bed_list[1]

    intersection = str((bed_1 + bed_2).count())
    b1_only = str((bed_1).count())
    b2_only = str((bed_2).count())
    
    print('\t'.join(["both", args.heads[0], args.heads[1]]))
    print('\t'.join([intersection, b1_only, b2_only]))


    generated_bed_intersection = bed_1 + bed_2 
    
    if args.o != None:
        remove_if_exists(args.o)
        with open(args.o, 'w+') as f:
            for line in generated_bed_intersection:
                f.write("\t".join(line))
                f.write("\n")
    else:
        pass

