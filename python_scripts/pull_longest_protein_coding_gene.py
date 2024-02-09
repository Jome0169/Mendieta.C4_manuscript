import argparse
import sys
import os
from Bio import SeqIO


"""
Parses protein file - and outputs the longest protein of a given isoform.
Assumes protein headers are in the form

>Zm00020ab330680_P002

Matches similar genes up utilizing the _P00 nomenclature, and prints out the
longest protein coding gene.

"""


def compare_dict_record_new_record(record, dict_record):
    if len(record.seq) > len(dict_record.seq):
        return True 
    elif len(record.seq) < len(dict_record.seq):
        return False 
    else:
        return False 

def read_fasta_file(input_file):

    gene_dict = {}

    for record in SeqIO.parse(input_file, "fasta"):
        gene = record.id.split("_")[0]

        if gene not in gene_dict:
            gene_dict[gene] = record
        elif gene in gene_dict:
            potential_replacment = compare_dict_record_new_record(record, gene_dict[gene])
            if potential_replacment == False:
                pass
            elif potential_replacment == True:
                gene_dict[gene] = record

    return gene_dict
        



def get_parser():
    parser = argparse.ArgumentParser(description='generated a sparse matrix for \
            later analysis')
    parser.add_argument('-input', '--input-file',  
            required=True, dest ="input")
    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()

    final_set = read_fasta_file(args.input)

    for key, val in final_set.items():
        print(">"+val.id)
        print(val.seq)
    







