"""
File: calculate_fpkm_over_region.py
Author: Pablo M
Email: yourname@email.com
Github: https://github.com/Jome0169
Description: Purpose of this script is to take in a HT-seq count file, as well
as the GTF file originally given to HT-seq count, and caluclate either the TPM
or for RPKM value. Note that FPKM at this moment is not supported as I have not
taken time to address fragment issues
"""
import argparse
import sys
import os
import itertools 


def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass


def read_in_bed(file_name):
    """TODO: Docstring for read_in_bed.
    :returns: TODO

    """
    line_list = []
    with open(file_name, 'r') as f:
        for line in f:
            clean_line = line.strip().split()
            line_list.append(clean_line)
    return line_list


def name_generator(base_name, list_len):
    """TODO: Docstring for name_generatory.
    :returns: TODO

    """
    counter = 1
    for item in range(0, list_len):
        updated_name = base_name + "_" + str(counter)
        yield updated_name 
        counter += 1

def generate_key_file(nested_list, base_name):
    """TODO: Docstring for generate_key_file.

    :nested_list: TODO
    :base_name: TODO
    :returns: TODO

    """
    key_file_dict = {}
    
    #Initiate Generator Function
    names = name_generator(base_name, len(nested_list))
    for item in nested_list:
        old_name = item[3]
        new_name = next(names)
        #print(old_name, new_name)
        #print([n for n in new_name])
        key_file_dict[old_name] = [old_name, new_name]
    return(key_file_dict)


def alter_name(loaded_bed_file, base_name, side):
    final_bed = []
    if side == "l":
        for item in loaded_bed_file:
            bed_name = item[3]
            new_bed_name = base_name  + "_" + bed_name
            item[3] = new_bed_name
            final_bed.append(item)

    elif side == "r":
        for item in loaded_bed_file:
            bed_name = item[3]
            new_bed_name = bed_name + "_" + base_name
            item[3] = new_bed_name
            final_bed.append(item)

    return(final_bed)

def write_key_file(value, file_name):
    """TODO: Docstring for write_key_file.
    :returns: TODO

    """
    with open(file_name, 'a+') as f:
        f.write('\t'.join(value))
        f.write('\n')



def write_output_bed(bed_list, output_file):
    """TODO: Docstring for write_key_file.
    :returns: TODO

    """
    if output_file != None:

        with open(output_file, 'a+') as f:
            for item in bed_list:
                f.write('\t'.join(item))
                f.write('\n')
    elif output_file == None:
        for item in bed_list:
            print('\t'.join(item))





def get_parser():
    parser = argparse.ArgumentParser(description='Calculates with the TPM or \
            FPKM given a gtf file, and a HT-seq count file output. ')
    parser.add_argument('-bed', "--bed_file", help="HT-seq Counts File to loadfile to load to \
            compute PKM values from", required=True, dest='bed')
    parser.add_argument('-base', "--gtf_file", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=False, dest='base')

    parser.add_argument('-m', "--mode", help="How to operate on bed name \
        present. Options are - replace - left - right \
        ", required=True, dest='mode')
    parser.add_argument('-key', "--key_file", help="Which value to \
            caluclate. OPtions: FPKM, TPM, RPKM", required = False, dest="key")
    parser.add_argument('-o', "--output", help="Bed file to load intervals. \
            These will be the intervals that will have their PKM value calcualte \
            ", required=False, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    if args.o != None:
        print("Working on %s" % args.bed)
        remove_if_exists(args.o)
    else:
        pass
    


    bed_file = read_in_bed(args.bed)

    if args.mode == "replace":
        remove_if_exists(args.key)
        key_file_dict = generate_key_file(bed_file, args.base)

        print("Working on %s" % args.bed)
        with open(args.o, 'a+') as f:
            for item in bed_file:
                #Write the key file
                write_key_file(key_file_dict[item[3]], args.key)
                #Grab the New referenced key, name new bed file
                item[3] = key_file_dict[item[3]][1]
                f.write('\t'.join(item))
                f.write('\n')

    elif args.mode == "left":
        modified_bed = alter_name(bed_file, args.base, 'l')
        write_output_bed(modified_bed, args.o)
    elif args.mode == "right":
        modified_bed = alter_name(bed_file, args.base, 'r')
        write_output_bed(modified_bed, args.o)
    else:
        sys.exit(-1)
