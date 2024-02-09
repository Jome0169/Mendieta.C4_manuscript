"""
Purpose of this script is to take in a blastp file where the entire input was a
series of genomic scaffolds, and the query database was a set of the longest
proteins in the maize genome. 

Input example (outfmt 6 from blast):

## qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
Zm00001eb093920_P002    chr7:185148527-185172422        47.059  85      45      0       22      106     5910    6164    1.74e-25        87.8
Zm00001eb093920_P002    chr7:185148527-185172422        77.778  27      6       0       101     127     6235    6315    1.74e-25        52.0
Zm00001eb093920_P002    chr9:124773081-124879240        42.748  131     48      1       22      125     454     846     4.91e-25        109

From this input file I need to weed out regions of the genome with signifigant
blast hits - take their location from the scaffold information provided here in
column 2 and "subtract" these regions - generating two new scaffolds. 

With all of this in mind the easiest fix is likely going to be using my old
friend pybedtools. 

First I'll read each coordinate region into a dictionary - with the list in the
value slot of all blast hits (which have passed filtering). From that I'll
itertivly go through each - generate a series of bedtools - and subtract them. 
"""

import argparse
import sys
import os




def read_blast_file(input_file):
    scaffold_with_hits_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split('\t')
            if cleaned_line[1] not in scaffold_with_hits_dict and float(cleaned_line[-2]) < .001:
                scaffold_with_hits_dict[cleaned_line[1]] = [cleaned_line]
            elif cleaned_line[1] in scaffold_with_hits_dict and float(cleaned_line[-2]) < .001:
                scaffold_with_hits_dict[cleaned_line[1]].append(cleaned_line)
            else:
                pass
    return(scaffold_with_hits_dict)


def parse_blast_dict(blast_dict):

    def generate_hit_region(scaffold_hit_string):
        """
        Takes the input string:
        chr7:185148527-185172422
        and returns [chr7, 185148527,185172422 ]
        """
        split_chrom_location = scaffold_hit_string.split(':')
        split_scaf_location = split_chrom_location[1].split("-")
        final_list = [split_chrom_location[0], int(split_scaf_location[0]), int(split_scaf_location[1])]
        return final_list

    def generate_updated_coordinates(chrom_location_1_location_2, blast_hit):
        """
        Takes the location [chr7,185148527,185172422] and adds the blast start
        and end locations from the blast_hist list to the start.

        This tells us the regoins of the original scaffold that needs to be
        removed.
        """
        blast_start = int(blast_hit[-4])
        blast_end = int(blast_hit[-3])


        if blast_start < blast_end:
            blast_scaffold_location_start = chrom_location_1_location_2[1] + blast_start
            blast_scaffold_location_end = chrom_location_1_location_2[1] + blast_end
        
        # Somtimes if the query hits on the negative strand - it flips the
        # scaffold coordinates. In which case switch the order 
        #['Zm00001eb308680_P001', 'chr7:79173296-79514962', '74.783', '115', '0', '1', '205', '290', '341599', '341255', '1.09e-33', '136']
        elif blast_start > blast_end:
            blast_scaffold_location_start = chrom_location_1_location_2[1] + blast_end
            blast_scaffold_location_end = chrom_location_1_location_2[1] + blast_start 
        
        final_return_bed_regions = [chrom_location_1_location_2[0],
                str(blast_scaffold_location_start),
                str(blast_scaffold_location_end)]

        return(final_return_bed_regions)
    
    finalized_dict = {}
    for region, blast_nested_list in blast_dict.items():
        finalized_dict[region] = []
        chrom_location_1_location_2 = generate_hit_region(region)
        ## qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        for blast_hit in blast_nested_list:
            blastp_bed_regoins = generate_updated_coordinates(chrom_location_1_location_2, blast_hit)
            finalized_dict[region].append(blastp_bed_regoins)
    return finalized_dict



def get_parser():
    parser = argparse.ArgumentParser(description='generated a sparse matrix for \
            later analysis')
    parser.add_argument('-i', '--input-file',  
            required=True, dest ="input")
    parser.add_argument('-o', '--output_base',
           required=False, dest ="output")
    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()
    
    blast_dict = read_blast_file(args.input)
    final_dict_to_subtract = parse_blast_dict(blast_dict)
    
    if args.output != None:
        with open(args.output, 'w+') as f:
            for key, val in final_dict_to_subtract.items():
                for region in val:
                    f.write('\t'.join(region))
                    f.write('\n')

    elif args.output == None:
        for key, val in final_dict_to_subtract.items():
            for region in val:
                print('\t'.join(region))




