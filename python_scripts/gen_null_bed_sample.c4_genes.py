import pybedtools
import argparse
import os
import itertools
import copy
import statistics
import csv
from threading import Thread
import random


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


def capture_gc_content(bed_file, genome_file):
    """Calls GC content using pybedtools. Requires the actual .fa genome.

    Returns the output from that 

    :arg1: TODO
    :returns: TODO

    """
    seq_profiled = bed_file.nucleotide_content(fi=genome_file).saveas()
    return(seq_profiled)


def calcualte_mean_GC_content(bed_file_nuc, start_index):

    GC_cont_all = [float(z[start_index + 1]) for z in bed_file_nuc]
    mean_gc_cont = statistics.mean(GC_cont_all)
    return(mean_gc_cont)


def calcualte_mean_counts(bed_file_nuc):
    """
    Calculates the mean count of intersections. The last integer here referes
    to the last column in the intersect command where this is called
    """

    mean_intersection = [int(z[-1]) for z in bed_file_nuc]
    calc_mean = statistics.mean(mean_intersection)
    return(calc_mean)


def calcualte_prop(bed_file_nuc):
    """
    Calculates the mean count of intersections. The last integer here referes
    to the last column in the intersect command where this is called
    """

    total_num = bed_file_nuc.count()
    
    final_num = bed_file_nuc.filter(lambda x: int(x[-1]) >= 1)
    count_final_num = final_num.count()

    proportion = count_final_num/total_num
    return(proportion)





def generate_matching_GC(acr_bed, count, genome_file, GC_mean, max_attempts=1000, current_attempt=1):
    """Shuffle features - and resample if features don't match mean GC content

    :bed_file: Bed file to shuffle
    :genome_file_index: You need this for shuffle, and unable to generate this
    within bedtools 

    """
    print(f"Attempt {current_attempt} out of {max_attempts}")

    if current_attempt > max_attempts:
        raise ValueError("Exceeded maximum attempts to sample matching GC content")

    #SHuffle these regions
    grab_equal_features = acr_bed.random_subset(count)
    
    #Count the number of fields present after this region survives (for use
    #later when calclating GC()
    surv_index = grab_equal_features.field_count()

    #Calc GC content
    rand_GC_content = capture_gc_content(grab_equal_features, genome_file)
    sample_mean_gc_score = calcualte_mean_GC_content(rand_GC_content,  surv_index)

    lower_bound = round(GC_mean * 0.95,3)
    upper_bound = round(GC_mean * 1.02,3)
    print(lower_bound, upper_bound)
    if lower_bound <= float(sample_mean_gc_score) <= upper_bound:
        return grab_equal_features
    else:
        return generate_matching_GC(acr_bed, count, genome_file, GC_mean, max_attempts, current_attempt + 1)
        
   

def capture_regions_with_matching_GC(acr_bed, number_matching_sites, genome_file, GC_mean):
    """Iteratively subsamples regions, checks their GC content and aggregates them if they match the desired GC content.

    :acr_bed: Bed file to subsample
    :genome_file: Required for pybedtools to get nucleotide content
    :GC_mean: Mean GC content to match against
    """
    # Calculate bounds for acceptable GC content
    lower_bound = round(GC_mean * 0.95, 3)
    upper_bound = round(GC_mean * 1.0, 3)

    # Initialize an empty pybedtools object to store matched regions
    matched_regions = []
    
    print("Capturing Matched Regions...")
    counter = 1
    while len(matched_regions) < number_matching_sites:
        print(f"Working....{counter}")
        print(f"Alread found {len(matched_regions)} regions, need {number_matching_sites - len(matched_regions)} more...")
        # Randomly subsample 1000 regions
        sampled_regions = acr_bed.random_subset(1000)
        sampled_regions_fields = sampled_regions.field_count()

        # Calculate the GC content for the subsampled regions
        regions_GC_content = capture_gc_content(sampled_regions, genome_file)

        for idx, region in enumerate(sampled_regions):
            region_gc_score = float(regions_GC_content[idx][sampled_regions_fields + 1])
            if lower_bound <= region_gc_score <= upper_bound and len(matched_regions) < number_matching_sites:
                #print("Found a passing region...")
                #print(idx, region, region_gc_score)
                matched_regions.append(region)
            else:
                pass
        counter += 1
    
    passing_region_bedtool = pybedtools.BedTool(matched_regions)
    return passing_region_bedtool




def keep_nfeatures(feature, n):
    """Given a list of features, keep only-
    n number of fields

    """
    new_feature = feature[:int(n)]
    feature = new_feature

    return feature




def filter_GC(feature, index_str, GC_num):
    GC_num_5= .10 * GC_num
    GC_num_high = GC_num + GC_num_5
    GC_num_low = GC_num - GC_num_5
    GC_score = float(feature[index_str + 1])

    if GC_score < GC_num_high and GC_score > GC_num_low:
        return(feature) 
    else:
        pass

def pull_not_passing_GC(feature, index_str, GC_num):
    GC_num_5= .10 * GC_num
    GC_num_high = GC_num + GC_num_5
    GC_num_low = GC_num - GC_num_5
    GC_score = float(feature[index_str + 1])
    
    if GC_score < GC_num_high and GC_score > GC_num_low:
        pass
    else:
        return(feature) 


def parse_gni(genome_index):
    chroms = {}
    with open(genome_index, 'r+') as f:
        for line in f:
            cleaned_line = line.strip().split()
            if cleaned_line[0] not in chroms:
                chroms[cleaned_line[0]] = (1,cleaned_line[1])
            else:
                pass
    return(chroms)


def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
        replicate peak calls, as well as unqiue peaks to each replicate and \
        outputs said peaks. ')
    parser.add_argument('-bed','--bed_file', help='Bed File to Mimic',\
        required=True, dest='bed'),
    parser.add_argument('-bed2','--bed_file2', help='Bed File to Mimic',\
        required=True, dest='bed2'),
    #parser.add_argument('-TFs','--TF_file', help='TF file in bed',\
    #    required=True, dest='TF_b'),
    #parser.add_argument('-exclusion_beds','--excl_files', help='Bed File to Mimic',\
    #    nargs="+", required=True, dest='exl'),
    parser.add_argument('-genome','--genome_file', help='Genome File to use', \
        required=True, dest='gn'),
    parser.add_argument('-genome_index','--genome_index', help='Genome File to use', \
        required=True, dest='gni'),
    parser.add_argument('-o','--output_name', help='output', 
        required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Read Files
    bed_file = read_in_bed_file(args.bed)
    ref_bed_file = read_in_bed_file(args.bed2)
    generate_genome_index_dic = parse_gni(args.gni)
    bed_file = bed_file.truncate_to_chrom(generate_genome_index_dic)
    ref_bed_file = ref_bed_file.truncate_to_chrom(generate_genome_index_dic)

    cell_type_split = {}
    ## Isoalte unique cell-types based off of the name in the 4th Col
    for i in ref_bed_file:
        grab_cell_type = i[3].split(';')[1]
        if grab_cell_type not in cell_type_split and "," not in grab_cell_type:
            cell_type_split[grab_cell_type] = [i]
        elif grab_cell_type in cell_type_split and "," not in grab_cell_type:
            cell_type_split[grab_cell_type].append(i)
        elif "," in grab_cell_type:
            pass
    
    ## Remove the broadly accessible class BC we'll be 
    ## using it as an outgroup of sorts
    grab_cell_types = list(cell_type_split.keys())
    grab_cell_types.remove("broadly_accessible")
    
    #All broadly acc acrs
    broadly_acc_acr_bed = pybedtools.BedTool(cell_type_split["broadly_accessible"]) 

    ## Calculate GC Content Per Cell Type
    gc_content_per_cell_type = {}
    bed_gc_content = capture_gc_content(bed_file, args.gn)
    take_len = bed_file.field_count()
    mean_cell_type_gc_score = calcualte_mean_GC_content(bed_gc_content, take_len)

    ## Generate matching list
    number_ACRs = bed_file.count()


    print("Taking equal sample of Broad ACRs...")
    matched_broad_acrs = capture_regions_with_matching_GC(broadly_acc_acr_bed, number_ACRs, args.gn, mean_cell_type_gc_score)


    #print("Taking equal sample of Broad ACRs...")
    #matched_broad_acrs = generate_matching_GC(broadly_acc_acr_bed,
    #                                          number_ACRs, args.gn, 
    #                                          mean_cell_type_gc_score, max_attempts = 2000)

    output_file_name = args.o + "." + ".broad_acr_null_list.bed"
    matched_broad_acrs.saveas(output_file_name)
