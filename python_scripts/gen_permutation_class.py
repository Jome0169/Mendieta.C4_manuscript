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





def generate_shuffled_features(new_bed, genome_file_index, genome_file, excl_regions, GC_mean):
    """Shuffle features - and resample if features don't match mean GC content

    :bed_file: Bed file to shuffle
    :genome_file_index: You need this for shuffle, and unable to generate this
    within bedtools 

    """
    

    rand_seed = random.randint(0,1000000)
    #SHuffle these regions
    shuffled_features = new_bed.shuffle(g = genome_file_index, excl = excl_regions,
            noOverlapping = True, seed = rand_seed).saveas()
    
    #Count the number of fields present after this region survives (for use
    #later when calclating GC()

    surv_index = shuffled_features.field_count()

    #Calc GC content
    rand_GC_content = capture_gc_content(shuffled_features, genome_file)
    

    #ID Proportion passing GC content Control - and failding
    not_passing_GC_content = rand_GC_content.each(pull_not_passing_GC, surv_index, GC_mean).saveas()
    passing_GC_content = rand_GC_content.each(filter_GC, surv_index, GC_mean).saveas()

    counter = 1
    while not_passing_GC_content.count() != 0:
    
        #Print message about how many are left
        print("Attempt Number %s to resample to pass GC" % str(counter))
        passing_GC_cont_count = passing_GC_content.count()
        not_passing_GC_cont_count = not_passing_GC_content.count()

        
        print("There are %s Passing GC Content Intitially" % str(passing_GC_cont_count))
        print("There are %s Not Passing GC Content Intitially" % str(not_passing_GC_cont_count))

        #Remove previous GC information
        remove_previous_GC = not_passing_GC_content.each(keep_nfeatures, 3).saveas()


        rand_seed = random.randint(0,1000000)
        reshuffle = remove_previous_GC.shuffle(g = genome_file_index, excl = excl_regions,
            noOverlapping = True, seed = rand_seed).saveas()

        
        #Reshuffle and count the files
        reshu_ind = remove_previous_GC.field_count()

        #Calc GC content
        reshuffle_GC = capture_gc_content(reshuffle, genome_file)

        #Filter into two groups
        recursive_not_passing_GC_content = reshuffle_GC.each(pull_not_passing_GC, reshu_ind, GC_mean).saveas()
        recursive_passing_GC_content = reshuffle_GC.each(filter_GC, reshu_ind, GC_mean).saveas()
        
        #Combine the passing regions together 
        passing_GC_content_2 = passing_GC_content.cat(recursive_passing_GC_content, force_truncate =
                True, postmerge=False).saveas()
        
        passing_GC_content = passing_GC_content_2
        passing_GC_content.saveas()

        not_passing_GC_content = recursive_not_passing_GC_content.saveas()
        counter += 1


    return(passing_GC_content)



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




def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
        replicate peak calls, as well as unqiue peaks to each replicate and \
        outputs said peaks. ')
    parser.add_argument('-bed','--bed_file', help='Bed File to Mimic',\
        required=True, dest='bed'),
    parser.add_argument('-TFs','--TF_file', help='TF file in bed',\
        required=True, dest='TF_b'),
    parser.add_argument('-exclusion_beds','--excl_files', help='Bed File to Mimic',\
        nargs="+", required=True, dest='exl'),
    parser.add_argument('-genome','--genome_file', help='Genome File to use', \
        required=True, dest='gn'),
    parser.add_argument('-genome_index','--genome_index', help='Genome File to use', \
        required=True, dest='gni'),
    parser.add_argument('-o','--output_name', help='output', \
        required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Read Files
    bed_file = read_in_bed_file(args.bed)
    
    #Regions Need to be greater than 10 bp for analysis
    bed_file_filtered = bed_file.filter(lambda x: len(x) >= 10)
    bed_file = bed_file_filtered.saveas()

    #Calc GC Content 
    gc_content = capture_gc_content(bed_file, args.gn)
    
    #Grab fields for index use
    take_len = bed_file.field_count()
    mean_gc_score = calcualte_mean_GC_content(gc_content, take_len)
    
    #We're going to have to filter these later for GC Content
    take_number_of_features = bed_file.count()
    

    #Read In Exclusion Regions. These are the regions which we should NOT put 
    exclusion_regions = [read_in_bed_file(i) for i in list(args.exl)]

    take_first_bed = exclusion_regions[0]
    rest_exlustion = exclusion_regions[1:]
    
    #Merge Exlusion Regions. Saved to a file for future use in shuffle which
    #only takes a file name and not a bedtool
    exclude_bed = args.o + ".exlude.bed"
    merged_exclusion_regions = take_first_bed.cat(*rest_exlustion, force_truncate = True, postmerge=False).sort().saveas(exclude_bed)
   
    
    #Intersect the real set of TFs with the readl set of CNSs. Save
    print("Intersecting Original set of TFs")
    base_o = args.o + ".intersecting_TF.counts.bed"
    target_region_TF_count = bed_file.intersect(args.TF_b, 
        c = True, sorted = True).saveas(base_o)            

    
    
    #We will write to this file as we go.
    generate_output_name = args.o + ".final.csv"  
    with open(generate_output_name, 'a+') as f:
        #f.write(",".join([str(z) for z in i]))
        #f.write('\n')

        #Calculate Mean Intersection of real intersect
        mean_true_intersection = calcualte_mean_counts(target_region_TF_count)
        prop_overlapping_TF = calcualte_prop(target_region_TF_count)
        print(prop_overlapping_TF)

        mean_intersection_vals = ["true_val", str(mean_true_intersection), str(prop_overlapping_TF)]

        f.write(",".join(mean_intersection_vals))
        f.write('\n')


        #Run the monte-carlo simulations 
        for i in range(1000):

            print("Iteration %s out of 1000" % str(i))
            #thread = Thread(target = generate_shuffled_features, args = (bed_file, args.gni, args.gn, exclude_bed, mean_gc_score))
            #thread.start()
            #shuffled_regions = thread.join().sort().saveas()
            
            print("Iteration %s out of 1000" % str(i))
            shuffled_regions = generate_shuffled_features(bed_file, args.gni, args.gn, exclude_bed, mean_gc_score).sort().saveas()
            
            print("Intersecting Null Intersection with TF File")
            null_intersection = shuffled_regions.intersect(args.TF_b, c = True, sorted = True)

            #number = random.randint(0,1000)
            #rand_file = "iter_file." + str(i) + ".bed"
            #null_intersection.saveas(rand_file)

            print("Take Mean of Null")
            mean_null_intersection = calcualte_mean_counts(null_intersection)
            prop_overlapping_TF_null = calcualte_prop(null_intersection)

            print("Appending to List")
            gen_strin = "iter_" + str(i)
            final_list = [gen_strin, str(mean_null_intersection), str(prop_overlapping_TF_null)]
            #mean_intersection_vals.append(final_list)


            f.write(",".join(final_list))
            f.write('\n')



   

    #average_len = statistics.mean([len(i) for i in bed_file])

    #x = pybedtools.BedTool()
    #generate_random_regions = x.random(l=average_len, g = args.gni).saveas()
    #
    #surviving_random_regions = generate_random_regions.intersect(merged_exclusion_regions, v=True, wa = True)
    #
    #surv_index = surviving_random_regions.field_count()
    #rand_GC_content = capture_gc_content(surviving_random_regions, args.gn)
    #
    #filtered_rand_regions = rand_GC_content.filter(filter_GC, surv_index, mean_gc_score).saveas()
    #random_subset_randomized_regions = filtered_rand_regions.random_subset(n=int(take_number_of_features)).sort().saveas()
    #
    #random_subset_randomized_regions_final = random_subset_randomized_regions.each(keep_nfeatures, 4).saveas()
    #

    #base_o = gs.o + ".intersecting_TF.counts.bed"

    #randomiz_o = args.o + ".randomization_regions.intersecting_TF.counts.bed"
    #
    #target_region_TF_count = bed_file.intersect(args.TF_b, 
    #        c = True, sorted = True).saveas(base_o)
    #filtered_rand_regions_TF_count = random_subset_randomized_regions_final.intersect(args.TF_b,
    #        c = True, sorted = True).saveas(randomiz_o)





