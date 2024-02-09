"""
The purpose of this script is to generate a genomic blacklist by either taking
in the INput ATAC-seq data, or Input ChIP-seq data. It accepts 4 column bed
files where the genome has been broken into 1kb bins, and the 4th column is
number of reads intersecting that region. 

If two files are given the mean of the bin is utilized.

This script after removing regions with low number of reads (zero inflated is
the norm), proceeds to take the reamining bins in the top 3rd quartile, and
reports regions within the top 1% of integration/enrichment. These regions
generally show substantial enzymatic bias, and are what we are attempting to
pull out.
"""

import argparse
import sys
import os
import pybedtools
import statistics
import numpy as np

def read_bed_file(bed_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    

    try:
        os.path.isfile(bed_file)
    except:
        FileNotFoundError

    bed_file_load = pybedtools.BedTool(bed_file).sort()    
    return bed_file_load

def replace_feature_name(feature, string_name):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    feature.name = string_name
    return feature

def trimmed_mean(data, percent):
    data = np.array(sorted(data))
    trim = int(percent*data.size/100.0)
    return data[trim:-trim].mean()


def take_median(bed_file, index):
    """
    Given a bed file, take the median of the score
    """

    filtered_bed = bed_file.filter(lambda x : float(x[index]) > 10)
    z = [float(i[index]) for i in filtered_bed]
    

    take_third_quartile_val = np.quantile(z, .75)
    
    #Only Take out things greater than third quartile
    gr_third_quartile = [x for x in z if float(x) > take_third_quartile_val]

    #Calculate median of third quartile
    bed_med = statistics.median(gr_third_quartile)

    top_1 = trimmed_mean(gr_third_quartile, 5)
    return(float(top_1))

def get_parser():
    parser = argparse.ArgumentParser(description='Calculate regions with fold \
            enrichment over X and outputs these regions.')
    parser.add_argument('-bed', "--bed_file", help="Column bed file with last \
    column being the scores to select on ", required=True, nargs = '+', dest='bed')
    parser.add_argument('-index', "--index_bed", help="Index to pull out", 
            required=True, dest='index')
    parser.add_argument('-fold', "--fold_dif", help=" Number of Times to fiter \
    over. 2 fold difference here would be 4. Must be integer ", \
    required=True, type = int, dest='fc')
    parser.add_argument('-o', "--output", help=" Output file. If not \
    given will be print ", required=False, dest='o')
    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    

    if len(args.bed) == 1:
        bed_file = read_bed_file(args.bed[0])

    elif len(args.bed) == 2:
        bed_file_1 = read_bed_file(args.bed[0])
        bed_file_2 = read_bed_file(args.bed[1])
        #bed_files_intersected = bed_file_1.intersect(bed_file_2, r=1)

        merged_beds = bed_file_1.cat(bed_file_2, postmerge=False).sort().saveas()

        bed_file = merged_beds.merge(c=4, o="mean", d=-1)

    elif len(args.bed) == 3:
        print("Cannot Work on Greater than 2 bed files. Exiting. Sorry :(")
        sys.exit(-1)


    median_all_bed = take_median(bed_file, int(args.index))
    generate_cut_off = float(args.fc) * median_all_bed

    passing_bed = bed_file.filter(lambda x: (float(x[int(args.index)])) > generate_cut_off).sort()




    if args.o != None:
        passing_bed.saveas(args.o)
    elif args.o == None:
        print(passing_bed)



