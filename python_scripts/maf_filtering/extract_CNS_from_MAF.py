"""
Purpose of this script is to take a MAF file with isolated CNS sequences and
output the orthologous regions in each of the species in the alignment. This
will use bipython as a method to deal with the data structure associated with
alignment.
"""


import argparse
import sys
import os
import pysam 
import pybedtools
import pandas as pd 

from Bio import pairwise2
from Bio import AlignIO
from Bio.Seq import Seq
import Bio.Align
#from Bio.AlignIO import MafIO



def isolate_conserved_sequences(multiple_alignment):
    keep_alignments = []
    for i in multiple_alignment:
        if len(i) > 1:
            alength = i.get_alignment_length()
            if alength > 5:
                keep_alignments.append(i)
    return(keep_alignments)


def calculate_prop_gap(sequence_algn):
    gaps = (sequence_algn.seq.count("-"))
    seq_len = (len(sequence_algn.seq))
    
    prop_missing = gaps/seq_len
    return(prop_missing)

     


def filter_alignment_on_qual(multiple_alignment):
    keep_alignments = []
    for i in multiple_alignment:
        #Grab the reference Sequences
        ref_seq_gap = calculate_prop_gap(i[0])
        print(i[0])
        print(ref_seq_gap)
        if ref_seq_gap > .5:
            print("Refrence sequence too gap filled - passing")
            print(i[0])
            break

        elif ref_seq_gap < .5:
            kept_sub_alignment = Bio.Align.MultipleSeqAlignment([])
            kept_sub_alignment.append(i[0])
            for sub_align in i[1:]:
                print(sub_align)
                align_prop_gap = calculate_prop_gap(sub_align)
                print(align_prop_gap)
                if align_prop_gap < .5:
                    #print(sub_align)
                    #print(align_prop_gap)
                    kept_sub_alignment.append(sub_align)
                elif align_prop_gap > .5:
                    pass
        keep_alignments.append(kept_sub_alignment)
    return(keep_alignments)

def generate_bed_line(sequence_algn, bed_CNS_id):
    """
    Given a MAF input - returns a complex list structure. First element is 
    species name, second element is a list of bed intervals for that given
    species, this facilitates iterative addition to bed dictionaries.
    """

    id_split = sequence_algn.id.split(".")
    species_name = id_split[0]
    chrom_name = id_split[-1]
    alignment_size = sequence_algn.annotations["size"]
    start = int(sequence_algn.annotations["start"]) + 1
    end = start + int(alignment_size)

    #bed_CNS_id =  "CNS_" + str(counter_val)
    bed_list_gen = [chrom_name, str(start), str(end), bed_CNS_id]
    return([species_name, bed_list_gen])


def generate_file_name(string):
    file_name = string + ".CNS_ID_regions.bed"
    return(file_name)

def generate_file_name_very_conserved(string):
    file_name = string + ".CNS_ID_regions.conserved_all_species.bed"
    return(file_name)





def remove_file_if_exits(file_name):
    try:
        os.remove(file_name)
    except OSError:
        pass



def get_parser():
    parser = argparse.ArgumentParser(description='Extract orthologous regions \
                                     associated with the MAF file for diverse species')
    parser.add_argument('-maf', "--maf_file", help="Maf file of alignment", 
                        required=True, dest='maf')
    parser.add_argument('-nspec', "--number_species", help="Number of \
                        alignments to expect for a CNS to be considered truly \
                        conserved", 
                    required=True, dest='nspec')
    #parser.add_argument('-fold', "--fold_dif", help=" Number of Times to fiter \
    #over. 2 fold difference here would be 4. Must be integer ", \
    #required=True, type = int, dest='fc')

    parser.add_argument('-odir', "--output_dir", help="Output directory. If not \
    given will write output files to current dir ", required=False, dest='o')
    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    
    #Parse MAF file
    maf_file = AlignIO.parse(args.maf, "maf")

    #Isoalte MAF alignments we want to keep
    isolated_cns = isolate_conserved_sequences(maf_file)

    print(len(isolated_cns))
    passing_CNS_sequences = filter_alignment_on_qual(isolated_cns)

    #Iterate and Annotate Alignments with CNS unique counters
    counter = 1
    species_bed_intersections = {}
    key_file = {}

    for align in passing_CNS_sequences:
        bed_CNS_id =  "CNS_" + str(counter)
        counter += 1

        if len(align) >= int(args.nspec):
            key_file[bed_CNS_id] = {}
            for sub_align in align:
                species_bed = generate_bed_line(sub_align, bed_CNS_id)
                sub_align.id = sub_align.id + ";" + bed_CNS_id

                if species_bed[0] not in species_bed_intersections:
                    species_bed_intersections[species_bed[0]] = [species_bed[1]]
                elif species_bed[0] in species_bed_intersections:
                    species_bed_intersections[species_bed[0]].append(species_bed[1])
                
                #Add the super conserved CNSs to a key file for referencing
                #later

                bed_string_safe = ".".join(species_bed[1])
                if species_bed[0] not in key_file[bed_CNS_id]:
                    key_file[bed_CNS_id][species_bed[0]] = bed_string_safe
                    #species_bed_intersections[species_bed[0]] = [species_bed[1]]
                elif species_bed[0] in key_file[bed_CNS_id]:
                    key_file[bed_CNS_id][species_bed[0]] = bed_string_safe
                    #species_bed_intersections[species_bed[0]].append(species_bed[1])
        else:
            for sub_align in align:
                species_bed = generate_bed_line(sub_align, bed_CNS_id)
                sub_align.id = sub_align.id + ";" + bed_CNS_id
                                                                                     
                if species_bed[0] not in species_bed_intersections:
                    species_bed_intersections[species_bed[0]] = [species_bed[1]]
                elif species_bed[0] in species_bed_intersections:
                    species_bed_intersections[species_bed[0]].append(species_bed[1])

        
        align.annotations["CNS_annot"] = bed_CNS_id
    
    remove_file_if_exits("Passing_CNS_alignments.maf")
    AlignIO.write(passing_CNS_sequences, "Passing_CNS_alignments.maf", "maf")


    key_file_very_conserved = pd.DataFrame.from_dict(key_file, orient = "index")
    key_file_very_conserved.to_csv('CNS_in_N_species.tsv', sep="\t")

    if args.o == None:
        for key in species_bed_intersections:
            species_x_name = generate_file_name(key)
            remove_file_if_exits(species_x_name)
            with open(species_x_name, 'a+') as f:
                for val in species_bed_intersections[key]:
                    all_strings = [str(i) for i in val]
                    f.write('\t'.join(all_strings))
                    f.write("\n")
            
            #Write only the stuff found in args.nspecies to different bed file
            species_name_very_conserved = generate_file_name_very_conserved(key)
            remove_file_if_exits(species_name_very_conserved)
            with open(species_name_very_conserved, 'a+') as f:
                for val in species_bed_intersections[key]:
                    all_strings = [str(i) for i in val]
                    if all_strings[3] in key_file:
                        f.write('\t'.join(all_strings))
                        f.write("\n")
                    else:
                        pass

    elif args.o != None:
        #Make Dir if does not exits
        if not os.path.exists(args.o):
            os.makedirs(args.o)
        else:
            pass
        
        #Write BED files per species
        for key in species_bed_intersections:
            species_x_name = generate_file_name(key)
            species_output_file_path = args.o + "/" + species_x_name
            remove_file_if_exits(species_output_file_path)
            with open(species_output_file_path, 'a+') as f:
                for val in species_bed_intersections[key]:
                    all_strings = [str(i) for i in val]
                    f.write('\t'.join(all_strings))
                    f.write("\n")

            #Write only the stuff found in args.nspecies to different bed file
            species_name_very_conserved = generate_file_name_very_conserved(key)
            species_very_conserved_output_file_path = args.o + "/" + species_name_very_conserved 
            remove_file_if_exits(species_very_conserved_output_file_path)
            with open(species_very_conserved_output_file_path, 'a+') as f:
                for val in species_bed_intersections[key]:
                    all_strings = [str(i) for i in val]
                    if all_strings[3] in key_file:
                        f.write('\t'.join(all_strings))
                        f.write("\n")
                    else:
                        pass





