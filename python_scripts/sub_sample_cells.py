"""
File: sub_sample_cells.py
Author: Pablo Mendieta 
Email: yourname@email.com
Github: https://github.com/yourname
Description: A quick script to subsample a bed file  based off a proportion of a meta
file. Takes N cells and outputs the integgration evens from all sampled cells
into a single bed file

Runs on multiple threads. Outputs bed file, as well as subsampled meta data. 
"""


import argparse
import sys
import os
import pybedtools
import pandas as pd
import random 
import string
import subprocess
from multiprocessing import Pool
from functools import partial

def read_bed_file(bed_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    print("Loading Bed file %s" % bed_file)
    

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

def split_bed_file(line, list_of_cell_IDs, generate_temp_file):
    """TODO: Docstring for split_bed_file.
    :returns: TODO

    """
    with open(generate_temp_file, 'a+') as f:
        if line.name in take_cell_names:
            f.write(str(line))
        else:
            pass


def get_parser():
    parser = argparse.ArgumentParser(description='Calculate regions with fold \
            enrichment over X and outputs these regions.')
    parser.add_argument('-bed', "--bed_file", help="Single Cell Bed file. cell \
            name needs to be in column name.", required=True, dest='bed'),
    parser.add_argument('-meta', "--meta_data", help="Meta data to read cells from", 
        required=True, dest='md'),
    parser.add_argument('-prop', "--proportion", help="Number of cells toa \
            sample", required=True, type = float, dest='prop'),
    parser.add_argument('-o', "--output", help=" Output file. If not \
        given will be print ", required=True, dest='o')
    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    print("Loading Meta Data")
    meta_data_read = pd.read_csv(args.md, sep='\t')

    if args.prop >= 1:
        print("Value greater than 1 - reverting to 1 to take 100% of cells")
        args.prop == 1
    else:
        pass

    #Subsample Rows
    subsampled_meta = meta_data_read.sample(frac=float(args.prop))
    take_cell_names = subsampled_meta["cellID"].tolist()

    count_number_reads = subsampled_meta["total"].sum()
    print(count_number_reads)

    print("Working on %s" % args.bed)
    bed_file = read_bed_file(args.bed)

    #Generate Random File Name for Storage
    #letters = string.ascii_lowercase
    #generate_temp_file = ''.join(random.choice(letters) for i in range(10))  + ".tmp"

    final_file_name = args.o + ".bed"
    print("Parsing Bed File. Please Hold.")
    with Pool(15) as pool:
        m = pool.map(
            partial(
                split_bed_file,
                list_of_cell_IDs = take_cell_names,
                generate_temp_file = final_file_name 
            ),
           bed_file)


    subsampe_meta_file_name= args.o + ".subsampled_meta.txt"
    subsampled_meta.to_csv(subsampe_meta_file_name, sep ='\t')
