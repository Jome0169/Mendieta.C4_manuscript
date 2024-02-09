"""
    new_peak_calling.call_ScACRs
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    This script takes single cell ATAC-seq data and call peaks on a per cluster
    basis. This script is an evolution of previous works in the schmitz lab
    done by Alex Marand and seeks to simplify this analyis by using a single
    scripts to run a series of proccesses and commands in order to ID cluster
    specific ACRs. 

    Methadology employed in this sciprt is taken from: 
    https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aav1898&file=aav1898_corces_sm.pdf
    
    And
    https://www.nature.com/articles/s41587-019-0206-z

    Example of Command: python call_ScACRs.py -bed testing_set.bed -meta Zm_leaf_svd.knn_100_strict.meta.annotation_V2.txt -col V1_annotation -base test -outdir TEST -fai Zm-B73-REFERENCE-NAM-5.0.chrom.size


    :copyright: (c) 2022 by Pablo Mendieta.
    :license: LICENSE_NAME, see LICENSE for more details.
"""

import argparse
import sys
import os
import pybedtools
import pandas as pd
import numpy
from multiprocessing import Pool, Manager
import multiprocessing
from functools import partial
import subprocess
import copy
import errno
import datetime
import random
import string


def read_bed_file(bed_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    try:
        os.path.isfile(bed_file)
        bed_file_load = pybedtools.BedTool(bed_file)
        return bed_file_load
    except:
        FileNotFoundError
        pass


def read_meta_file(meta_file):
    """Read meta data in using pandas. Tab needs to seeperate meta file
    :returns: TODO

    """
    read_meta = pd.read_csv(meta_file, sep="\t")
    return read_meta


def remove_file(file_name):
    """General purpose remove file function called in downstrem functions
    :returns: NA
    """

    try:
        os.remove(file_name)
    except OSError:
        pass


def replace_feature_name(feature, string_name):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    feature.name = string_name
    return feature


def replace_feature_name_w_coords(feature):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    new_feature = feature
    genenrate_feature_name = (
        str(feature.chrom)
        + "__"
        + str(feature.start)
        + "__"
        + str(feature.stop)
        + "__"
        + str(feature.score)
    )
    new_feature.name = genenrate_feature_name
    #feature.name = genenrate_feature_name
    return new_feature


def mkdir_p(path):
    """
    If output path is given attempts to write this output directory to store
    files
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def extend_summits(feature):
    """Bedtools pydfunction. Calculated per kb accessability.
    :returns: Appends per kb accessability as the last index of the string.

    """
    new_start = feature.start - 250
    new_stop = feature.stop + 250
    feature.start = new_start
    feature.stop = new_stop
                                                                            
    if new_start < 0:
        pass
    else:
        return feature


def regen_interval_name(list_1):
    """Regenerate bed interval with correct score after splitting.

    :arg1: TODO
    :returns: TODO

    """

    generate_string = list_1[0] + "__" +  list_1[1] + "__" +  list_1[2] + "__" +  list_1[3]

    final_list = list_1
    final_list[3] = generate_string
    return(final_list)

def rip_grep_split(i, bed_file_name, rep_1_dict, rep_2_dict, cell_file_names):
    """Splits bed file by cell type into pseudo replicates. This is critical as
    it allows us to test for signifigance further downstream.

    :arg1: TODO
    :returns: TODO

    """
    def write_quick_temp(list_1, temp_file_name):
        """TODO: Docstring for write_quick_temp.
        :returns: TODO

        """
        with open(temp_file_name, 'a+') as f:
            for i in list_1:
                f.write(i)
                f.write('\n')

        
    #Random letters here are to allow threading
    letters = string.ascii_lowercase
    rand_string = ''.join(random.choice(letters) for i in range(3)) 
    tmp_file_1 = i + "_" + rand_string + "_sampled_cells_rep1.txt"
    tmp_file_2 = i + "_" + rand_string + "_sampled_cells_rep2.txt"
    tmp_pool_file = i + "_" + rand_string +"_sampled_cells_pool.txt"

    remove_file(tmp_file_1)
    remove_file(tmp_file_2)
    remove_file(tmp_pool_file)

    write_quick_temp(rep_1_dict[i], tmp_file_1)
    write_quick_temp(rep_2_dict[i], tmp_file_2)

    gathered_cell_list = rep_1_dict[i] + rep_2_dict[i]
    write_quick_temp(gathered_cell_list, tmp_pool_file)

    output_file_rep_1 = cell_file_names[i][1]
    output_file_rep_2 = cell_file_names[i][2]
    output_file_pool = cell_file_names[i][0]

    generate_rg_command_3 = f"grep -Ff '{tmp_pool_file}' {bed_file_name} > {output_file_pool}"
    generate_rg_command_1 = f"grep -Ff '{tmp_file_1}' {output_file_pool} > {output_file_rep_1}"
    generate_rg_command_2 = f"grep -Ff '{tmp_file_2}' {output_file_pool} > {output_file_rep_2}"

    print(f"Running grep Command {generate_rg_command_3}")
    subprocess.run([generate_rg_command_3], shell=True)
    print(f"Running grep Command {generate_rg_command_1}")
    subprocess.run([generate_rg_command_1], shell=True)
    print(f"Running grep Command {generate_rg_command_2}")
    subprocess.run([generate_rg_command_2], shell=True)

    remove_file(tmp_file_1)
    remove_file(tmp_file_2)
    remove_file(tmp_pool_file)






def simple_split(bed_file, cellID_celltype_dict, rep_1_dict, rep_2_dict):
    """TODO: Docstring for simple_split.
    :returns: TODO

    """
    cell_ID_dict_info = dict()
    for i in (bed_file):
        cell_name = list(i)[3]
        if cell_name in cellID_celltype_dict:
            grab_cell_type = cellID_celltype_dict[cell_name]
            if grab_cell_type not in cell_ID_dict_info:
                cell_ID_dict_info[grab_cell_type] = [[], []]
            else:
                pass
            if cell_name in rep_1_dict:
                cell_ID_dict_info[grab_cell_type][0].append(i)
            elif cell_name in rep_2_dict:
                cell_ID_dict_info[grab_cell_type][1].append(i)
        else:
            pass
    return cell_ID_dict_info



def write_pooled_bed(meta_data, bed_file, column_name, output_base_name, replicate_col, output_dir, cores):
    """TODO: Docstring for write_pooled_bed.

    :meta_data: TODO
    :cell_types: TODO
    :output_base_name: TODO
    :output_dir: TODO
    :: TODO
    :returns: TODO


    Generated a dictionary with cell_type being the key, and the various output
    files being a list. For example: 
    ['testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.pool.bed', 'testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.rep1.bed', 'testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.rep2.bed', 'testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.pool.macs', 'testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.rep1.macs', 'testing/epidermis_ncell_3607/Zm.peaks.annot_v4_no_ground_epidermis_ncell_3607.rep2.macs']

    This list is iteravivley added to and used throughout the script in order
    to reference and generate various outputs files for use downstream.

    """

    def generate_cellID_celltype_dict(pandas_df, column_name):
        """Geneartes a dictionary where key is cellID from meta data and value
        is the annotation label. 
        :returns: dictionary where the cell ID is the key and the annotation is
        the value

        """
        cellID_celltype_dict = {}
        for i in pandas_df[["cellID", column_name]].itertuples():
            id_name = i[1]
            cell_annot = i[2]
            cellID_celltype_dict[id_name] = cell_annot
        return cellID_celltype_dict



    def generate_celltype_cellID_dict(pandas_df, column_name):
        """Geneartes a dictionary where key is annotation_lable from meta data and value
        is the cell ID. 
        :returns: dictionary where the cell ID is the key and the annotation is
        the value
                                                                                
        """
        celltype_cellID_dict = {}
        for i in pandas_df[["cellID", column_name]].itertuples():
            id_name = i[1]
            cell_annot = i[2]
            if cell_annot not in celltype_cellID_dict:
                celltype_cellID_dict[cell_annot] = [id_name]
            elif cell_annot in celltype_cellID_dict:
                celltype_cellID_dict[cell_annot].append(id_name)
        return celltype_cellID_dict


    def generate_file_names(list_of_annots, output_base_name, output_dir):
        """Generates file output names to be used downstream for the file
        split. Will generate an output nested list 
        :returns: cell_annotation_1 : [annotation_type_1.pool, annotation_type_1.rep1,
        annotation_type_1.rep2], etc...]

        """
        file_name_nested_list = {}

        for i in list_of_annots:
            if output_dir == None:
                mkdir_p(i)
                pool = i + "/" + output_base_name + "_" + i + ".pool.bed"
                rep1 = i + "/" + output_base_name + "_" + i + ".rep1.bed"
                rep2 = i + "/" + output_base_name + "_" + i + ".rep2.bed"
                file_name_nested_list[i] = [pool, rep1, rep2]
            elif output_dir != None:
                mkdir_p(output_dir + "/" + i)
                pool = (
                    output_dir
                    + "/"
                    + i
                    + "/"
                    + output_base_name
                    + "_"
                    + i
                    + ".pool.bed"
                )
                rep1 = (
                    output_dir
                    + "/"
                    + i
                    + "/"
                    + output_base_name
                    + "_"
                    + i
                    + ".rep1.bed"
                )
                rep2 = (
                    output_dir
                    + "/"
                    + i
                    + "/"
                    + output_base_name
                    + "_"
                    + i
                    + ".rep2.bed"
                )
                file_name_nested_list[i] = [pool, rep1, rep2]
        return file_name_nested_list

    def write_output_bed_pools_reps(
        cell_type_interval_dictionary, cell_type_file_names
    ):
        """TODO: Docstring for write_output_bed_pools_reps.
        :returns: TODO

        """
        for cell_type, interval_pools in cell_type_interval_dictionary.items():

            ### Write the replicates
            print(f"Writing output bed files for {cell_type}")
            grab_pool_name = cell_type_file_names[cell_type][0]
            grab_rep_1_name = cell_type_file_names[cell_type][1]
            grab_rep_2_name = cell_type_file_names[cell_type][2]

            pybedtools.BedTool(interval_pools[0]).sort().saveas(grab_rep_1_name)
            pybedtools.BedTool(interval_pools[1]).sort().saveas(grab_rep_2_name)

            pooled_intervals = interval_pools[0] + interval_pools[1]

            pybedtools.BedTool(pooled_intervals).sort().saveas(grab_pool_name)

    cell_types = pd.unique(meta_file[column_name])

    ## THe Below Sectoin generates a series of dictionaries with cellIDs as
    ## keys and the val field being the actual associated cell type.
    ## These are eventually used with the dictionary of cellType ; File name

    ##Generate Data For InterGroup Splits
    split_by_cell_type = meta_file.groupby(by=[column_name])
    size_values = meta_file.groupby(by=[column_name]).size()

    cellID_annotation_type_rep_1_dictionary = {}
    cellID_annotation_type_rep_2_dictionary = {}

    annotation_type_cellID_rep_1_dictionary = {}
    annotation_type_cellID_rep_2_dictionary = {}
    if replicate_col == None:
        for i in cell_types:
            ##Grab and Shuffle
            grabbed_cellID = meta_file.loc[meta_file[column_name] == i].sample(frac=1)
            group_len = round(len(grabbed_cellID) / 2)

            rep_1_group = grabbed_cellID.head(group_len)
            rep_2_group = grabbed_cellID.tail(group_len)

            dict_sub_group_rep_1 = generate_cellID_celltype_dict(rep_1_group, column_name)
            dict_sub_group_rep_2 = generate_cellID_celltype_dict(rep_2_group, column_name)

            grouping_val_dict_sub_group_rep_1 = generate_celltype_cellID_dict(rep_1_group, column_name)
            grouping_val_dict_sub_group_rep_2 = generate_celltype_cellID_dict(rep_2_group, column_name)

            cellID_annotation_type_rep_1_dictionary.update(dict_sub_group_rep_1)
            cellID_annotation_type_rep_2_dictionary.update(dict_sub_group_rep_2)

            annotation_type_cellID_rep_1_dictionary.update(grouping_val_dict_sub_group_rep_1) 
            annotation_type_cellID_rep_2_dictionary.update(grouping_val_dict_sub_group_rep_2) 

    elif replicate_col != None:
        replicate_values = meta_file[replicate_col].unique().tolist()
        print(f"Using the following replicates to splot {replicate_values}")

        if len(replicate_values) > 2:
            print(f"Expected two replicates, reverting to psuedo replicate \
                  methods. Sorry, improve the code on your own :) ")
            sys.exit(-1)
        elif len(replicate_values) == 2:
            for i in cell_types:
                grabbed_cellID = meta_file.loc[meta_file[column_name] == i].sample(frac=1)

                rep_1_group = grabbed_cellID[grabbed_cellID[replicate_col] == replicate_values[0]]

                rep_2_group = grabbed_cellID[grabbed_cellID[replicate_col] == replicate_values[1]]


                print(rep_1_group.head())
                print(rep_2_group.head())
                                                                                                            
                dict_sub_group_rep_1 = generate_cellID_celltype_dict(rep_1_group, column_name)
                dict_sub_group_rep_2 = generate_cellID_celltype_dict(rep_2_group, column_name)
                                                                                                            
                grouping_val_dict_sub_group_rep_1 = generate_celltype_cellID_dict(rep_1_group, column_name)
                grouping_val_dict_sub_group_rep_2 = generate_celltype_cellID_dict(rep_2_group, column_name)
                                                                                                            
                cellID_annotation_type_rep_1_dictionary.update(dict_sub_group_rep_1)
                cellID_annotation_type_rep_2_dictionary.update(dict_sub_group_rep_2)
                                                                                                            
                annotation_type_cellID_rep_1_dictionary.update(grouping_val_dict_sub_group_rep_1) 
                annotation_type_cellID_rep_2_dictionary.update(grouping_val_dict_sub_group_rep_2) 

    cellID_celltype_dict = generate_cellID_celltype_dict(meta_data, column_name)
    cell_type_file_name = generate_file_names(cell_types, output_base_name, output_dir)

    x = [remove_file(i) for i in cell_type_file_name]

    print(f"Preparing to Split Bed file By Cell Type Using {cores} Cores")
    gather_cell_types_to_split = annotation_type_cellID_rep_1_dictionary.keys()
    with Pool(cores) as pool:
        m = pool.map(
            partial(
                rip_grep_split,
                bed_file_name = bed_file,
                rep_1_dict = annotation_type_cellID_rep_1_dictionary,
                rep_2_dict = annotation_type_cellID_rep_2_dictionary,
                cell_file_names=cell_type_file_name,
            ),
           gather_cell_types_to_split 
        )


    
    return cell_type_file_name


def sub_func_macs2(i, gsize, dict_file_name, output_dir):
    """Using the pooled and split bed files go ahead and call peaks using
    MACS2. Note that the qvalue here is set and further corrected downstream.
    
    :arg2: 
    :returns: TODO
    
    """
    # dict_file_names_macs2_output = copy.deepcopy(dict_file_names)
    # for cell_file_name in dict_file_name[i]:

    # cell_peak_calls = dict_file_name[i][-1]
    # cell_file_name = dict_file_name[i][1]

    grab_file_list = dict_file_name[i]
    length = len(grab_file_list)
    middle_index = length // 2

    read_bed_files = grab_file_list[:middle_index]

    ## This is kinda nasty, but takes the output MACS2 file name which has the
    ## path appended, splits on the "/" and then returns the actual file
    output_file_names = [x.split("/")[-1] for x in grab_file_list[middle_index:]]

    ##Take the first three files, and following three files associaated with the
    ##output file named rep.MACS
    if output_dir == None:
        final_output_dir_name = i
    elif output_dir != None:
        final_output_dir_name = output_dir + "/" + i

    for bed_f, out_f in zip(read_bed_files, output_file_names):
        # cell_peak_calls = cell_file_name[-1]
        if output_dir == None:
            generate_macs2_command = f"macs2 callpeak -t {bed_f} -f BED -g {gsize} --nomodel \
            --keep-dup auto --extsize 150 --shift -50 --qvalue .05 --outdir {final_output_dir_name} --bdg \
            -n {out_f}"
        elif output_dir != None:
            generate_macs2_command = f"macs2 callpeak -t {bed_f} -f BED -g {gsize} --nomodel \
            --keep-dup auto --extsize 150 --shift -50 --qvalue .05 --outdir {final_output_dir_name} --bdg \
            -n {out_f}"

        print(f"Running MACS2 Command {generate_macs2_command}")
        subprocess.run([generate_macs2_command], shell=True, check=True)

        print("Done Running MACS2 Calls")



def run_macs2_threaded(dict_file_name, gsize, output_directory, cores):
    """Calls the function to run MACS2 on both pool. This allows for rapid
    processing of multiple cell-types using a threaded option since MACs2 is
    pretty efficient
    :returns: TODO

    """

    def generate_macs_file_name(dict_file_names, output_dir):
        dict_file_names_macs2_output = copy.deepcopy(dict_file_names)
        # Deep Copy Dit
        for cell_type, cell_file_names in dict_file_names.items():
            for cell_file_name in cell_file_names:
                cell_peak_calls = cell_file_name.split("/")[-1].replace(".bed", ".macs")

                if output_dir == None:
                    final_output_dir_name = cell_type
                elif output_dir != None:
                    final_output_dir_name = output_dir + "/" + cell_type

                save_output_file = final_output_dir_name + "/" + cell_peak_calls
                dict_file_names_macs2_output[cell_type].append(save_output_file)

        return dict_file_names_macs2_output

    updated_dict_file_names = generate_macs_file_name(dict_file_name, output_directory)
    grab_cell_keys = updated_dict_file_names.keys()
    for key, val in updated_dict_file_names.items():
        print(key, val)

    with Pool(int(cores)) as pool:
        m = pool.map(
            partial(
                sub_func_macs2,
                gsize = gsize,
                dict_file_name=updated_dict_file_names,
                output_dir=output_directory,
            ),
            grab_cell_keys,
        )

    return updated_dict_file_names



def read_reference_genome_file(fai_file):
    """TODO: Docstring for read_reference_genome_file.
    :returns: TODO

    """
    chrom_size = {}
    with open(fai_file, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split()
            if cleaned_line[0] not in chrom_size:
                chrom_size[cleaned_line[0]] = int(cleaned_line[1])
            else:
                pass
    return(chrom_size)



def filter_bedgraph(line, chrom_size_dict):
    """TODO: Docstring for fitler_bedgraph.
    :returns: TODO

    """
    grab_chrom = line[0]
    if int(line[2]) < chrom_size_dict[grab_chrom]:
        return(line)
    else:
        return(None)

def generate_normalized_bw(file_dict, reference_file, output_dir):
    """TODO: Docstring for generate_normalized_bw.

    :arg1: TODO
    :returns: TODO

    """

    chrom_size_dict = read_reference_genome_file(reference_file)
    
    updated_dict_file_names = copy.deepcopy(file_dict)
    for cell_type, file_list in file_dict.items():
        bed_pool_file = file_list[0]
        bed_graph_file = file_list[3] + "_treat_pileup.bdg"
        bed_graph_file_sorted = file_list[3] + "_treat_pileup.sorted.bdg"
        bed_graph_file_normalized = file_list[3] + "_treat_pileup.normalized.bdg"
        bw_file_name = file_list[3].replace(".pool.macs", ".normalized.bw")

        bed_pool_file_read = read_bed_file(bed_pool_file)
        number_reads_in_pool = bed_pool_file_read.count()
        generate_CPM_value = number_reads_in_pool / 1000000

        print(
            f"The Scaling value for {cell_type} is {generate_CPM_value} based off of {number_reads_in_pool}"
        )

        try:
            subprocess.run(
                [f"sort -k 1,1 -k2,2n {bed_graph_file} > {bed_graph_file_sorted}"],
                shell=True,
                check=True,
            )
        except ValueError:
            print("Something in the sort Failed")
            exit(-1)

        normalized_CPM_values = []
        with open(bed_graph_file_sorted, "r") as f:
            for line in f:
                cleaned_line = line.strip().split("\t")
                normalized_value = round(
                    (float(cleaned_line[-1]) / generate_CPM_value), 5
                )
                cleaned_line[-1] = str(normalized_value)
                normalized_CPM_values.append(cleaned_line)

        remove_file(bed_graph_file_normalized)
        ##Write potential removal if exits
        with open(bed_graph_file_normalized, "a+") as f:
            for line in normalized_CPM_values:
                filtered_line = filter_bedgraph(line,chrom_size_dict)
                if filtered_line != None:
                    f.write("\t".join(line))
                    f.write("\n")
                else:
                    pass

        try:
            generate_bedgraph_to_bw = f"bedGraphToBigWig {bed_graph_file_normalized} {reference_file} {bw_file_name}"
            print(f"Running BedGraphToBigWig Command {generate_bedgraph_to_bw}")
            subprocess.run([generate_bedgraph_to_bw], shell=True)
        except subprocess.CalledProcessError as e:
            print(
                "Error in Running Bedgraph to Bigwig. Is UCSC Loaded? Files corrupted. \
                    Conntinuing on in the pipeline."
            )



def filter_peaks_for_fragments(peak_file, tn5_file):
    """
    Remove peaks which appear to be caused by single fragments
    """

    save_peak_file = peak_file.saveas()
    #Count the original number of fields and isolate where the coverage tract is
    #going to be 
    coverage_slot = save_peak_file.field_count() + 1
    
    #calcualte the coverage of the given peaks and filter based off of number
    #of Tn5 
    coverage_calculations = save_peak_file.coverage(tn5_file).saveas()

    filtered_peaks = coverage_calculations.filter(
        lambda x: int(x[coverage_slot]) >= 20 
    ).saveas()

    return(filtered_peaks)


def take_reproducible_peaks_shuffle(file_dict, index_fai_file, output_dir):
    """TODO: Docstring for take_reproducible_peaks.
    :returns: TODO

    """
    
    #def coverage_fileter(feature, slot):
    #    updated_feauture = 


    updated_dict_file_names = copy.deepcopy(file_dict)

    for cell_type, file_list in file_dict.items():
        tn5_file = read_bed_file(file_list[0]).sort()

        #Generate file names for reference and use
        pool_peak_file = read_bed_file(file_list[3] + "_peaks.narrowPeak")
        rep1_peak_file = read_bed_file(file_list[4] + "_peaks.narrowPeak")
        rep2_peak_file = read_bed_file(file_list[5] + "_peaks.narrowPeak")

        reproducible_peaks = file_list[3].replace(
            ".pool.macs", ".reproducible_narrow_peaks"
        )
        null_shuffle_peaks = file_list[3].replace(".pool.macs", ".null_shuffle")
        repro_integration = file_list[3].replace(".pool.macs", ".integration.tn5")
        null_shuffle_peaks_inter = file_list[3].replace(
            ".pool.macs", ".null.integration.tn5"
        )
        reproducible_summits = file_list[3].replace(
            ".pool.macs", ".reproducible_summits"
        )
    
        print(f"Filtering Peaks for {cell_type} by Removing Those Generated by a single Fragment")
        #Filtering summits for fragment issues
        pooled_summit_peaks = read_bed_file(file_list[3] + "_summits.bed")
        pooled_summit_peaks_extended = pooled_summit_peaks.saveas().each(extend_summits)
        pool_peaks_filtered = filter_peaks_for_fragments(pooled_summit_peaks_extended, tn5_file)


        print(f"Generating Reproducible Peaks for {cell_type}")
        reproducible_peaks_local = (
            pool_peak_file.intersect(pool_peaks_filtered, u=True, f=0.1)
            .intersect(rep1_peak_file, u=True, f=0.1)
            .intersect(rep2_peak_file, u=True, f=0.1)
            .sort()
            .saveas(reproducible_peaks)
        )

        print(f"Generating null dist of Peaks for {cell_type}")
        generated_shuffle = (
            reproducible_peaks_local.shuffle(excl=reproducible_peaks, g=index_fai_file)
            .sort()
            .saveas(null_shuffle_peaks)
        )

        print(f"Generating Reproucible Summits {cell_type}")
        pooled_summit_peaks.intersect(reproducible_peaks_local, u=True).sort().saveas(reproducible_summits)

        print(f"Interseting pooled bed with  null dist and repro Peaks for {cell_type}")
        intersection_real_peaks = reproducible_peaks_local.intersect(
            tn5_file, c=True, sorted=True
        ).saveas(repro_integration)
        intersection_fake_peaks = generated_shuffle.intersect(
            tn5_file, c=True, sorted=True
        ).saveas(null_shuffle_peaks_inter)

        updated_dict_file_names[cell_type].append(reproducible_peaks)
        updated_dict_file_names[cell_type].append(null_shuffle_peaks)
        updated_dict_file_names[cell_type].append(repro_integration)
        updated_dict_file_names[cell_type].append(null_shuffle_peaks_inter)
        updated_dict_file_names[cell_type].append(reproducible_summits)
    return updated_dict_file_names


def calculate_per_kb_accessability(feature):
    """Bedtools pydfunction. Calculated per kb accessability.
    :returns: Appends per kb accessability as the last index of the string.

    """
    calculated_vale = int(feature[-2]) / (
        (int(feature.stop) - int(feature.start)) / 1000
    )
    # feature[-1] = str(calculated_vale)

    new_feature = feature
    new_feature[-1] = str(calculated_vale)
    return new_feature


def extend_fields(feature, n):
    fields = feature.fields[:]
    while len(fields) < n:
        fields.append(".")
    return pybedtools.create_interval_from_list(fields)


def eFDR_process(dict_file_names, fdr_val, output_dir):
    """Set null distribution of peaks based off of null peak integrations. From
    here set FDR and call peaks which pass threshold.
    :returns: TODO

    """

    updated_dict_file_names = copy.deepcopy(dict_file_names)
    for key, val in dict_file_names.items():
        reproducible_peaks = read_bed_file(val[-3])
        null_values = read_bed_file(val[-2])

        #try:
        final_field_count_repro = reproducible_peaks.field_count() + 1
        final_field_count_null = null_values.field_count() + 1
        repro_tn5_denisty = reproducible_peaks.each(
            extend_fields,final_field_count_repro 
        ).each(calculate_per_kb_accessability).saveas()
        null_values_density = null_values.each(
            extend_fields, final_field_count_null 
        ).each(calculate_per_kb_accessability).saveas()
        
        try:
            take_null_tn5_quantile = numpy.quantile(
                [float(x[-1]) for x in null_values_density], 1 - fdr_val
            )
            print(f"This Is the null Tn5 Quantile {take_null_tn5_quantile}")
        except:
            take_null_tn5_quantile = 10

        final_output_name = val[-3].replace(".integration.tn5", ".passing.integration.tn5")

        filtered_real_sites = repro_tn5_denisty.filter(
            lambda x: float(x[-1]) > take_null_tn5_quantile
        ).saveas(final_output_name)
        peaks_passing = filtered_real_sites.count()
        print(f"There are {peaks_passing} passing filtering steps")

        updated_dict_file_names[key].append(final_output_name)

    return updated_dict_file_names


def remove_overlapping_peaks(py_bed_tool):
    """TODO: Docstring for remove_overlapping_peaks.
    :returns: TODO

    """
    
    

    def ID_correct_interval(arg1):
        """TODO: Docstring for ID_correct_interval.

        :arg1: TODO
        :returns: TODO

        """
        split_to_lists = arg1[3].split(",")
        split_intervals_to_string = [i.split("__") for i in split_to_lists]
        sorted_lists = sorted(
            split_intervals_to_string, key=lambda x: max(x[-1]), reverse=True
        )
        most_signifigant = sorted_lists[0]
        sig_intersect = regen_interval_name(most_signifigant)
        generated_interval = pybedtools.create_interval_from_list(sig_intersect)
        return generated_interval

    final_list = []
    for i in py_bed_tool:
        if "," in i[3]:
            most_sig_interval = ID_correct_interval(i)
            final_list.append(most_sig_interval)
        elif "," not in i[3]:
            final_list.append(i)

    generated_bed_file = pybedtools.BedTool(final_list)
    return generated_bed_file


def grab_summits_normalize_beds(dict_file_names, output_dir):
    """TODO: Docstring for grab_summits_normalize_beds.
    :returns: TODO

    """

    def generate_normalized_acc_score(bed_file):
        total = 0
        for i in bed_file:
            total += float(i[-1])
        final_value = total / 1000000
        return final_value

    def normalize_values(feature, norm_val):
        """TODO: Docstring for normalize_values.
    
        :arg1: TODO
        :returns: TODO
    
        """
        take_freature_score = float(feature[-1])
        normalized_score = take_freature_score / norm_val
        feature[-1] = str(normalized_score)
        return feature


    


    updated_dict_names = copy.deepcopy(dict_file_names)
    for key, val in dict_file_names.items():
        repro_summit_name = val[3].replace(
            ".pool.macs", ".reproducible_summits.passing_FDR"
        )

        read_repro_peak_file = read_bed_file(val[-2])
        read_summit_file = read_bed_file(
            val[3].replace(".pool.macs", ".reproducible_summits")
        )

        intersecting_summits = read_summit_file.intersect(
            read_repro_peak_file, u=True
        ).saveas()

        # Gather Acc score per million
        per_peak_norm_value = generate_normalized_acc_score(intersecting_summits)

        # Update score and extend summits +/- 250
        final_norm_summits = (
            intersecting_summits.each(normalize_values, per_peak_norm_value)
            .each(extend_summits)
            .saveas()
            .each(remove_invaid_merge)
            .saveas()
            .each(replace_feature_name_w_coords)
            .saveas()
            .each(remove_invaid_merge)
            .saveas()
        )

        try:
            # Fitler the overlapping peaks to the most signifigant peaks
            merged_signigifigant_peaks = final_norm_summits.merge(
                c=4, o="collapse"
            ).each(remove_invaid_merge).saveas()
            most_signifigant_peaks = remove_overlapping_peaks(
                merged_signigifigant_peaks
            )
            most_signifigant_peaks.each(remove_invaid_merge).saveas(repro_summit_name)
        except:
            print(f"{key} failed - moving on ")
            with open(repro_summit_name, "a+") as f:
                pass

        updated_dict_names[key].append(repro_summit_name)

    return updated_dict_names


def remove_invaid_merge(feature):
    """ Given a bedtool feature, go through and replace the feature name. All
    Features will have the same name.
    """
    if feature.start < feature.stop and feature.start > 0:
        return feature
    elif feature.start < feature.stop and feature.start < 0:
        pass
    elif feature.start > feature.stop:
        pass
    else:
        pass


def add_final_counter(bed_tool, base_name):
    """TODO: Docstring for add_final_counter.
    :returns: TODO

    """
    final_bed_tool = []
    counter = 1
    for i in bed_tool:
        new = i
        ## Grab the Normalized Score
        grabbed_score = i[3].split("_")[-1]
        generate_new_name = base_name + "_" + str(counter)
        new[3] = generate_new_name
        new[4] = grabbed_score
        final_bed_tool.append(new)
        counter += 1
    return final_bed_tool

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
    parser = argparse.ArgumentParser(
        description="Call Peaks for scATAC data. \
    Requires cluster annnotations, as well as BED file ipput."
    )
    parser.add_argument(
        "-bed",
        "--bed_file",
        help="4 Column bed file with last \
    column being the scores to select on ",
        required=True,
        dest="bed",
    )
    parser.add_argument(
        "-meta",
        "--meta_data",
        help="Meta_data to grab things \
    from",
        required=True,
        dest="m",
    )
    parser.add_argument(
        "-col",
        "--column_annotation",
        help="Column to call \
    peaks on from",
        required=True,
        dest="col",
    )
    parser.add_argument(
        "-rep",
        "--replicate_col",
        help="Replicate column to split on. If not included psuedo reps \
        are generated for your cells",
        required=False,
        dest="rep",
    )
    parser.add_argument(
        "-bw",
        "--bigwig",
        help="Generate CPM normalized bigwig file per cluster",
        required=False,
        dest="bw",
    )
    parser.add_argument(
        "-fdr",
        "--fdr",
        help="Set the FDR rate for peak callling. Default is .001",
        required=False,
        dest="fdr",
    )

    parser.add_argument(
        "-cores",
        "--cores",
        help="Number of cores to run the script with",
        required=False,
        dest="cores",
    )

    parser.add_argument(
        "-gsize",
        "--genome_size",
        help="Genome Size to pass to macs2",
        required=False,
        dest="gsize",
    )
    parser.add_argument(
        "-fai", "--fai_file", help="Chrom size file", required=False, dest="fai",
    )

    parser.add_argument(
        "-base", "--base_name", help="Output basename.", required=True, dest="o"
    )
    parser.add_argument(
        "-outdir", "--output_dir", help="Output directory.", required=False, dest="od"
    )

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":


    args = get_parser().parse_args()
    meta_file = read_meta_file(args.m)
    print("Working on %s" % args.bed)
    bed_file = read_bed_file(args.bed)


    generate_genome_index_dic = parse_gni(args.fai)

    #Setting the number of cores for splitting
    if args.cores == None:
        cores = 5
    elif args.cores != None:
        try:
            cores = int(args.cores)
        except:
            print(f"Your number of cores {args.cores} cannot be converted to a \
                  int. Canceling")
            exit(1)

    if args.rep == None:
        print(f"No replicate column given, cells will be psuedo rep'd and ACRs \
              will be called. Note - here be technical dragons")
        args.rep = None
    if args.rep != None:
        pass



    ## Split bed file and write output
    pooled_bed_file_dict = write_pooled_bed(meta_file, args.bed, args.col, args.o, args.rep, args.od, cores)
    
    if args.gsize == None:
        print(f"Setting genome size to maize default 1.6e9 for peak calling.... You \
              will want to fix this.... ")
    elif args.gsize != None:
        pass

    ## Call Peaks using MACS2
    called_macs2_peaks = run_macs2_threaded(pooled_bed_file_dict, args.gsize, args.od, cores)
    
    ## Generating Bigwigs if Requested. Requires and FAI file to be given if
    ## flag given
    if args.bw == None:
        pass
    elif args.bw != None and args.fai == None:
        print(
            f"Unable to generate BW files, a fai file must be provided if \
        generating BWs. Cancelling BW Creation. \n Continuing"
        )
        pass
    elif args.bw != None and args.fai != None:
        pass
        generate_normalized_bw(called_macs2_peaks, args.fai, args.od)
    

    # Generate a Null distribution to generate FDR comparisons 
    reproducible_null = take_reproducible_peaks_shuffle(
        called_macs2_peaks, args.fai, args.od
    )
    
    
    if args.fdr == None:
        fdr_rate = .001
    elif args.fdr != None:
        try:
            fdr_rate = float(args.fdr)
        except: 
            print(f"Your FDR rate of {args.fdr} cannot be converted to a float. Canceling")
            exit(1)

    ## Running FDR Analysis
    threhold_passed_beds = eFDR_process(reproducible_null, fdr_rate, args.od)

    ## Remove overlapping peaks within each cell type. Take most signifigant
    final_dict_list = grab_summits_normalize_beds(threhold_passed_beds, args.od)
    
    ## Gather all passing peaks per cell type and add them to the same bed file
    all_passing_500_bp_peaks = [x[-1] for x in final_dict_list.values()]
    read_all_500_bp_peaks = [read_bed_file(z) for z in all_passing_500_bp_peaks]

    if len(read_all_500_bp_peaks) > 1:
        cated_bed_files = (
            read_all_500_bp_peaks[0]
            .cat(*read_all_500_bp_peaks[1:], postmerge=False, force_truncate=False)
            .saveas()
        )
    elif len(read_all_500_bp_peaks) == 1:
        cated_bed_files = read_all_500_bp_peaks[0].saveas()


    ## Remove potential BS issues associated with non-real overlaps
    cated_bed_files_cleaned = cated_bed_files.each(remove_invaid_merge).sort().saveas()

    ## Merge peaks from different cell  types select most signifigant
    merged_cleaned_bed_files = cated_bed_files_cleaned.merge(c=4, o="collapse").each(remove_invaid_merge).truncate_to_chrom(generate_genome_index_dic).saveas()

    print(f"Finishing Up Peak Identification")
    generate_final_peak_set = remove_overlapping_peaks(merged_cleaned_bed_files)
    numbered_acr_peaks = pybedtools.BedTool(generate_final_peak_set).sort()
    number_acr_fields = numbered_acr_peaks.field_count() + 1
    numbered_acr_peaks_extended = numbered_acr_peaks.each(
        extend_fields, number_acr_fields
    )

    ## Adding counter to final ACR names
    numbered_acr_peaks_final = add_final_counter(numbered_acr_peaks_extended, "scACR")

    if args.od == None:
        output_file_name = args.o + ".500bp_peaks.bed"
    elif args.od != None:
        output_file_name = args.od + "/" + args.o + ".500bp_peaks.bed"

    print(f"Output peaks are in {output_file_name}")
    pybedtools.BedTool(numbered_acr_peaks_final).truncate_to_chrom(generate_genome_index_dic).saveas(output_file_name)
