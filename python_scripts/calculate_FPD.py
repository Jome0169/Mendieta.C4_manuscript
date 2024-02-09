import pyBigWig
import argparse
import sys
import os
import copy 
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

    bed_file_load = pybedtools.BedTool(bed_file)
    return bed_file_load

def read_bw(bw_file):
    """Reads in the big wig file provided using hte pybigwig librayr from
    deepTools

    :bw_file: TODO
    :returns: TODO

    """
    bw = pyBigWig.open(bw_file)

    if bw.isBigWig() == True:
        return(bw)
    elif bw.isBigWig() == False:
        sys.exit("Failed Bigwig")


def add_feat_name(feature, string_name):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    new_name = feature.name + string_name

    feature.name = new_name
    return feature

def extend_fields(feature, n):
    fields = feature.fields[:]
    while len(fields) < n:
        fields.append('.')
    return pybedtools.create_interval_from_list(fields)

def add_name_bed(bed_file):
    """

    :bed_file: TODO
    :returns: TODO

    """
    take_len = range(1,len(bed_file))
    final_list = []

    for read, number in zip(bed_file, take_len):
        new_name = read[3] + '_number_' + str(number)
        read[3] = new_name

        final_list.append(read)


    final_bed_tool = pybedtools.BedTool(final_list)
    return(final_bed_tool)




def FDP_regions(feature, n, direction, chrom_dict):
    """Given a list of features, expand them by N bp up or downstream.
    """

    grab_chrom = feature.chrom
    chrom_len = chrom_dict[grab_chrom]

    if direction == "up":
        new_start = int(feature.start) - n
        if new_start < 0:
            new_start = 0
        else:
            pass

        new_stop = int(feature.start) 

    elif direction == "down":
        new_start = int(feature.stop) 
        new_stop = int(feature.stop) + n

        if new_stop > chrom_len:
            new_stop = chrom_len
        else:
            pass
    #Alter Feature
    feature.start = new_start
    feature.stop = new_stop
    return feature

def cal_mean_bw(feature, bw):
    """Given a list of features, expand them by N bp up or downstream.
    """
    chrom = feature.chrom
    start = int(feature.start)
    stop = int(feature.stop) 
    
    scored_mean = bw.stats(str(chrom), start, stop, type="mean")

    if scored_mean != None:
        feature[-1] = str(scored_mean[0])
    elif scored_mean == None or len(scored_mean) == 0 :
        feature[-1] = str(0)

    return feature




def calculate_FPD(motif_bed, upstream_bed, downstream_bed):
    """Calculates the FPD values using the bed tools

    :motif_bed: TODO
    :upstream_bed: TODO
    :downstream_bed: TODO
    :returns: TODO

    """
    score_dictionary ={}

    for read in motif_bed:
        if read[3] not in score_dictionary:
            score_dictionary[read[3]] = [read[-1]]
        elif read[3] in score_dictionary:
            print("ERROR MOTIF BEING CALCULATED TWICE")

    for upstream_read in upstream_bed:
        if upstream_read[3] in score_dictionary:
            score_dictionary[upstream_read[3]].append(upstream_read[-1])
        else:
            print("ERROR READ NOT THERE")
            pass

    for downstream_read in downstream_bed:
        if downstream_read[3] in score_dictionary:
            score_dictionary[downstream_read[3]].append(downstream_read[-1])
        else:
            pass
    return(score_dictionary)


def calculate(score_dict):
    """TODO: Docstring for calculate.
    :returns: TODO

    """
    final_dict = copy.deepcopy(score_dict)

    for key,val in final_dict.items():

        try:
            flank_1 = float(val[1])
        except ValueError:
            flank_1 = 0 
        
        try:
            flank_2 = float(val[2])
        except ValueError:
            flank_2  = 0

        try:
            motif = float(val[0])
        except ValueError:
            motif = 0

        mean_flanks = np.mean([flank_1, flank_2])

        FPD = motif - mean_flanks

        final_dict[key].append(str(FPD))
    return(final_dict)



def read_genome_file(file_name):
    """TODO: Docstring for read_genome_file.
    :returns: TODO

    """
    final = {}
    with(open(file_name, 'r')) as f:
        for line in f:
            cleaned_line = line.strip().split()
            if cleaned_line[0] not in final:
                final[cleaned_line[0]] = int(cleaned_line[1])
            else:
                pass
    return(final)


def get_parser():
    parser = argparse.ArgumentParser(description='Calculate regions with fold \
            enrichment over X and outputs these regions.')
    parser.add_argument('-bed', "--bed_file", help="Column bed file with last \
    column being the scores to select on ", required=True, dest='bed')
    parser.add_argument('-bw', "--bw", help="Bw with integration events", 
            required=True, dest='bw')
    parser.add_argument('-g', "--genome_file", help=" Number of Times to fiter \
    over. 2 fold difference here would be 4. Must be integer ", \
    required=True, dest='g')

    parser.add_argument('-s', "--size", help=" Number of Times to fiter \
    over. 2 fold difference here would be 4. Must be integer ", \
    required=True, dest='s')
    parser.add_argument('-o', "--output", help=" Output file. If not \
    given will be print ", required=False, dest='o')
    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    genome_file_read = read_genome_file(args.g) 
    #Set up bed file
    bed_file_raw = read_bed_file(args.bed)
    bed_file_renamed = add_name_bed(bed_file_raw)
    final_field_count = bed_file_renamed.field_count() + 1

    bed_file = bed_file_renamed.each(extend_fields, final_field_count).remove_invalid().saveas()
    

    #Read big wig
    bw = read_bw(args.bw)
    
    #Generate upstreamd and downstream flanks
    upstream_regions = bed_file.each(FDP_regions, int(args.s),
            "up",genome_file_read).each(extend_fields, final_field_count).remove_invalid().saveas()

    downstream_regions = bed_file.each(FDP_regions, int(args.s),
            "down",genome_file_read).each(extend_fields, final_field_count).remove_invalid().saveas()

    
    #Score regions
    scored_motif = bed_file.each(cal_mean_bw, bw).remove_invalid().saveas()
    scored_upstream = upstream_regions.each(cal_mean_bw, bw).remove_invalid().saveas()
    scored_downstream = downstream_regions.each(cal_mean_bw, bw).remove_invalid().saveas()
    
    #Calcualte FPD
    dict_of_scores = calculate_FPD(scored_motif, scored_upstream, scored_downstream)
    calculate_FPD_score = calculate(dict_of_scores)

    
    #Output
    if args.o != None:
        pass

    elif args.o == None:
        headers = ["region", "motif_score", "upstream_flank",
                "downstream_flank", "FPD"]

        joined_headers = '\t'.join(headers)
        print(joined_headers)
        for key,val in calculate_FPD_score.items():

            key_list = [key]
            final_list = key_list + val
            joined_vals = '\t'.join(final_list)
            print(joined_vals)


