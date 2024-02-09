from sys import argv
import argparse
import sys
import os

def remove_file(file_name):
    """TODO: Docstring for remove_file.
    :returns: TODO

    """

    try:
        os.remove(file_name)
    except OSError:
        pass


collapse_dict = { "columella_meristem" : "columella",
"inflorescence_meristem": "L1",
"L1_branch_meristem" : "L1",
"L1_glume_primordia" : "L1",
"L1_inflorescence_meristem" : "L1",
"L1_leaf_primordia" : "L1",
"L1_lemma_primordia" : "L1",
"L1_lower_floral_meristem" : "L1",
"L1_palea_primordia" : "L1",
"L1_pistil_primordia" : "L1",
"L1_proximal_meristem" : "L1",
"L1_root_tip" : "L1",
"L1_SAM" : "L1",
"L1_spikelet_meristem" : "L1",
"L1_spikelet_pair_meristem" : "L1",
"L1_stamen_primordia" : "L1",
"L1_upper_floral_meristem" : "L1",
"protophloem_sieve_element": "phloem",
"metaphloem_sieve_element": "phloem",
"phloem_pole_pericycle": "phloem",
"phloem_sieve_element_precursors": "phloem",
"pith_parenchyma" : "parenchyma",
"phloem_parenchyma": "parenchyma",
"parenchyma" : "parenchyma",
"metaxylem" : "xylem",
"protoxylem" : "xylem",
"axillary_meristem": "meristem",
"procambial_meristem": "meristem",
"proximal_meristem": "meristem",
"branch_meristem": "meristem",
"ground_meristem": "meristem",
"lower_floral_meristem": "meristem",
"upper_floral_meristem": "meristem",
"root_apical_meristem": "meristem",
"spikelet_meristem": "meristem",
"spikelet_pair_meristem": "meristem",
"leaf_primordia" : "primordia",
"glume_primordia" : "primordia",
"lemma_primordia" : "primordia",
"palea_primordia" : "primordia",
"secondary_root_primordia" : "primordia"}

def write_split_markers(lines_to_split, tis_type, marker_base_file):
    """Takes all lines from marker file and splits commans write_split_markers.
    if multiple commas - splits on them, and writes each line/tisuse sepretely

    """
    with open(marker_base_file, 'a') as f:
        for item in lines_to_split:
            comm_count = item[-1].count(',')
            if comm_count >= 1:
                take_base_list = item[0:5]
                split_comma = item[-1].split(',')
                for cell_type in split_comma:
                    final_list = take_base_list + [cell_type]
                    catch_bool = purge_lines(final_list, tis_type, )
                    if catch_bool != None:
                        f.write('\t'.join(catch_bool))
                        f.write('\n')
                    else:
                        pass

            elif comm_count == 0:
                catch_bool = purge_lines(item, tis_type, )
                if catch_bool != None:
                    f.write('\t'.join(catch_bool))
                    f.write('\n')
                else:
                    pass



def generate_compressed_marker_list(lines_to_split, collapse_dict):
    """TODO: Docstring for generate_compressed_marker_list.
    :returns: TODO

    """
    gathered_collapsed_markers = []
    for item in lines_to_split:
        comm_count = item[-1].count(',')

        if comm_count >= 1:
            take_base_list = item[0:5]
            split_comma = item[-1].split(',')
            for cell_type in split_comma:
                if cell_type in collapse_dict:
                    grab_new_name = collapse_dict[cell_type]
                    final_list = take_base_list + [grab_new_name]
                    gathered_collapsed_markers.append(final_list)
                elif cell_type not in collapse_dict:
                    final_list = take_base_list + [cell_type]
                    gathered_collapsed_markers.append(final_list)

        elif comm_count == 0:
            cell_type = item[-1]
            take_base_list = item[0:5]
            if cell_type in collapse_dict:
                grab_new_name = collapse_dict[cell_type]
                final_list = take_base_list + [grab_new_name]
                gathered_collapsed_markers.append(final_list)
            elif cell_type not in collapse_dict:
                final_list = take_base_list + [cell_type]
                gathered_collapsed_markers.append(final_list)

    unique_collapsed_only = []

    for item in gathered_collapsed_markers:
        if item not in unique_collapsed_only:
            unique_collapsed_only.append(item)

    return(unique_collapsed_only)

def write_compressed_markers(lines_to_write, tis_type, compressed_base_name):
    """Writes  the compressed versionn of the marker file .write_compressed_markers.
    :Written file using the above compressed_base_name dict

    """
    with open(compressed_base_name, 'a') as f:
        for item in lines_to_write:
            catch_bool = purge_lines(item, tis_type)
            if catch_bool != None:
                f.write('\t'.join(catch_bool))
                f.write('\n')
            else:
                pass


def purge_lines(line_to_process, tissue_type):
    """TODO: Docstring for purge_lines.

    :line_to_process: TODO
    :exlusion_dict: TODO
    :returns: TODO

    """
    exlusion_dict = {"leaf": ["inflorescence", "root", "pistil", "spikelet",
        "atrichoblast", "trichoblast", "locule", "glume", "exodermis",
        "lemma", "floral", "branch_meristem", "axillary_meristem"
        "lodicule", "stamen", "endodermis", "pericycle"],

        "ear": ["root", "endodermis", "pericycle","pistil", "stamen", "leaf"],

        "root": ["inflorescence", "pistil", "spikelet", "stamen", "leaf", "bundle_sheath", 
            "developing_pavement_cell",'L1']}

    grab_correct_exclusion_dict = exlusion_dict[tissue_type]
    take_base_name = line_to_process[-1]

    if all(s not in take_base_name for s in grab_correct_exclusion_dict):
        return(line_to_process)
    else:
        return(None)



def get_parser():
    parser = argparse.ArgumentParser(description='Given a list of markers \
    prepares  them to be used in the R script annotate_cells.v3.R.')

    parser.add_argument('-base', "--base_file",  help="base name to be used in the script", 
            required=True, dest='base')
    parser.add_argument('-tis', "--tis_name",  help="Dict used to exlude specific markers and terms", 
        required=True, dest='tis')
    parser.add_argument('-marker', "--marker", \
            help="File of markers to split",  \
            required=True, dest='mark')

    parser.add_argument('-header', "--header", \
            help="Add a header line? Yes? No?",  \
            required=False, action='store_true', dest='hd')



    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    
    base_name = args.base + ".split_markers.txt"
    base_name_compressed = args.base + ".compressed_markers.txt"
    
    remove_file(base_name)
    remove_file(base_name_compressed)
    
    if args.hd == True:
        header_line = ["chr","start","end","geneID", "name", "type"]
        cleaned_file_lines = [header_line]
    elif args.hd == False:
        cleaned_file_lines = []

    with open(args.mark, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split("\t")
            cleaned_file_lines.append(cleaned_line)



    write_split_markers(cleaned_file_lines, args.tis, base_name)
    collapsed_markers = generate_compressed_marker_list(cleaned_file_lines, collapse_dict)
    write_compressed_markers(collapsed_markers, args.tis, base_name_compressed)




