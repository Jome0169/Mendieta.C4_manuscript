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
        if len(i) >= 2:
            alength = i.get_alignment_length()
            if alength > 8:
                keep_alignments.append(i)
    return(keep_alignments)


def calculate_prop_gap(sequence_algn):
    gaps = (sequence_algn.seq.count("-"))
    seq_len = (len(sequence_algn.seq))
    prop_missing = gaps/seq_len
    return(prop_missing)

     


def filter_alignment_on_qual(multiple_alignment):
    keep_alignments = []
    len_alingmnet = len(multiple_alignment)

    counter = 0
    for i in multiple_alignment:
        #Grab the reference Sequences
        ref_seq_gap = calculate_prop_gap(i[0])
        if ref_seq_gap > .5:
            counter += 1
            pass
        elif ref_seq_gap < .5:
            kept_sub_alignment = Bio.Align.MultipleSeqAlignment([])
            kept_sub_alignment.append(i[0])
            for sub_align in i[1:]:
                align_prop_gap = calculate_prop_gap(sub_align)
                if align_prop_gap < .5:
                    kept_sub_alignment.append(sub_align)
                elif align_prop_gap > .5:
                    pass

            keep_alignments.append(kept_sub_alignment)

    print(f"In total {counter} alignemnts were discarded due to too many gaps \
          out of a total of {len_alingmnet} alignmtns")
    return(keep_alignments)

def generate_bed_line(sequence_algn):
    """
    Given a MAF input - returns a complex list structure. First element is 
    species name, second element is a list of bed intervals for that given
    species, this facilitates iterative addition to bed dictionaries.

    Sbicolor-Sorghum_bicolor.Chr05(+)/65900685-65900727      
    Pmiliaceum-Pmiliaceum.CM009707.2(-)/5687297-5687339      
    Ufusca-Ufusca.Chr08(+)/22490049-22490091                 
    Osativa-Osativa.Chr11(+)/20138252-20138294
    """

    first_split = sequence_algn.id.split("/")
    
    grab_first_half = first_split[0]
    if "." in grab_first_half:
        split_on_dot = grab_first_half.split(".")
        if len(split_on_dot) == 2:
            chrom_name = split_on_dot[1].split("(")[0]
            species_name = split_on_dot[0].split("-")[0]
        elif len(split_on_dot) > 2:
            chrom_name = ".".join(split_on_dot[1:3]).split("(")[0]
            species_name = split_on_dot[0].split("-")[0]
    elif "." not in grab_first_half:
        split_hypen = grab_first_half.split("-")
        chrom_name = split_hypen[-1].split("(")[0]
        species_name = "-".join(split_hypen[0:2])


    split_coords = first_split[1].split(("-"))
    
    start = split_coords[0]
    split_last_segment = split_coords[1].split(';')
    
    
    
    end = split_last_segment[0]
    bed_CNS_id = split_last_segment[1]
    
    bed_list_gen = [chrom_name, str(start), str(end), bed_CNS_id]
    return([species_name, bed_list_gen])

def update_alignment_name(sub_alignment, CNS_ID_string):
    updated_sub_id = sub_alignment.id + ";" + CNS_ID_string
    sub_alignment.id = updated_sub_id
    return sub_alignment

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
    maf_file = AlignIO.parse(args.maf, "clustal")

    #Isoalte MAF alignments we want to keep
    isolated_cns = isolate_conserved_sequences(maf_file)
    remove_low_quality_alignments = filter_alignment_on_qual(isolated_cns)
    passing_CNS_sequences = isolate_conserved_sequences(remove_low_quality_alignments)

    # This list is used later to populate the dictionary in case a conserved
    # CNS doesn't have all species. Prevents key filr from being messed up.
    all_species = []

    
    #For each alignment, rename the alignment with a unique CNS ID value

    key_file = {}
    counter = 1
    for align in passing_CNS_sequences:
        bed_CNS_id =  "CNS_" + str(counter)
        if len(align) >= int(args.nspec):
            #These are the really conserved groups of CNSs
            #Add them to the key dictionary for later processing
            key_file[bed_CNS_id] = {}
            for sub_align in align:
                rename_ID = sub_align.id + ";" + bed_CNS_id
                sub_align.id = rename_ID
            counter += 1
        elif len(align) < int(args.nspec) and len(align) > 0:
            for sub_align in align:
                rename_ID = sub_align.id + ";" + bed_CNS_id
                sub_align.id = rename_ID
            counter += 1
        else:
            pass


    species_bed_intersections = {}
    for align in passing_CNS_sequences:
        for sub_align in align:
            species_bed = generate_bed_line(sub_align)

            #Populate list to fix later
            if species_bed[0] not in all_species:
                all_species.append(species_bed[0])
            
            #Add the CNS to the species bed intersection no matter what
            if species_bed[0] not in species_bed_intersections:
                species_bed_intersections[species_bed[0]] = [species_bed[1]]
            elif species_bed[0] in species_bed_intersections:
                species_bed_intersections[species_bed[0]].append(species_bed[1])
            
            #If CNS in the key file of ultra conserved stuff, add a safe
            #version of the bed intersection to the dictionary
            grab_bed_CNS = species_bed[1][-1]
            if grab_bed_CNS in key_file:
                bed_string_safe = ".".join(species_bed[1])
                if species_bed[0] not in key_file[grab_bed_CNS]:
                    key_file[grab_bed_CNS][species_bed[0]] = bed_string_safe
                elif species_bed[0] in key_file[grab_bed_CNS]:
                    key_file[grab_bed_CNS][species_bed[0]] = bed_string_safe
            else:
                pass


    remove_file_if_exits("Passing_CNS_alignments.maf")
    remove_file_if_exits("CNS_in_N_species.tsv")
    AlignIO.write(passing_CNS_sequences, "Passing_CNS_alignments.maf", "maf")

    
    for key, val in key_file.items():
        for species_name in all_species:
            if species_name not in val:
                val[species_name] = "NA"
            elif species_name in val:
                pass

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





