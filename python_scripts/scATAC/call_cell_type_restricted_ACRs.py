from sys import argv
import argparse
import sys
import os
import pandas as pd
import numpy as np 
import pybedtools
from tqdm import tqdm

def remove_file(arg1):
    """TODO: Docstring for remove_file(.

    :arg1: TODO
    :returns: TODO

    """
    try:
        os.remove(arg1)
    except OSError:
        pass

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


def read_DA_file(da_file_name, remove_start_name, remove_tail_name, filter_value):
    """TODO: Docstring for read_DA_file.
    :returns: TODO

    """
    
    

    if remove_start_name != None and remove_tail_name != None:
        updated_key_name = da_file_name.replace(remove_start_name, "").replace(remove_tail_name, "")
        if "ncell" in updated_key_name:
            base_key_name = updated_key_name.split("_ncell_")[0]
        elif "ncell" not in updated_key_name:
            base_key_name = updated_key_name
        else:
            base_key_name = updated_key_name
    else:
            base_key_name = da_file_name

    
    #Read in the file and label which cell type it came from
    DA_read_file = pd.read_csv(da_file_name, sep='\t', lineterminator='\n')
    passing_DA_regions = DA_read_file[DA_read_file["padj"] < float(filter_value)]
    passing_DA_regions["cell_type"] = base_key_name
    
    #return the final name to call the cell types, as well as the filtered DA 
    #List of genes
    z = [base_key_name, passing_DA_regions]
    return(z)


def combine_ct_ACRs(df_list):
    """TODO: Docstring for list_of_DA_ACR_dfs.
    :returns: TODO

    """
    combined_ACR_DA_values = pd.concat(df_list)
    sub_sampled_df = combined_ACR_DA_values[["gene_name", "cell_type"]]
    #unique_acrs = sub_sampled_df[["gene_name"]].drop_duplicates()
    #print(unique_acrs)

    #final = pd.merge(unique_acrs, sub_sampled_df, on="gene_name", how = "left")
    merged_cell_types_ACRs = sub_sampled_df.groupby(["gene_name"])["cell_type"].apply(list).reset_index()
    #merged_cell_types_ACRs.columns = ["gene_name", "cell_type"]
    #merged_cell_types_ACRs.set_index("gene_name")
    return(merged_cell_types_ACRs)



def write_output(list_n, output_name):
    """TODO: Docstring for write_output.
    :returns: TODO

    """
    remove_file(output_name)
    with open(output_name, 'a+') as f:
        for line in list_n:
            f.write(str(line))




def get_parser():
    parser = argparse.ArgumentParser(description='Given a list of markers \
    prepares  them to be used in the R script annotate_cells.v3.R.')
    parser.add_argument('-DA', "--DA_files",  help="Differential Accessability Files", 
            required=True, dest='DA', nargs ="+")

    parser.add_argument('-bed', "--bed_file",  help="Bed file of ACRs", 
        required=True, dest='bed')
    parser.add_argument('-pval', "--pvalue",  help="Adjusted P value to filter on", 
        required=True, dest='pval', type = float)

    parser.add_argument('-front', "--tfront",  help="Adjusted P value to filter on", 
        required=False, dest='fr', type = str)
    parser.add_argument('-tail', "--ttail",  help="Adjusted P value to filter on", 
        required=False, dest='tl', type = str)

    parser.add_argument('-o', "--output", help="Output File",  \
            required=False, dest='o')

    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    

    #read in all DA 
    DA_cell_types = {}
    for i in args.DA:
        x = read_DA_file(i, args.fr, args.tl, args.pval)
        DA_cell_types[x[0]] = x[1]

    
    all_gathered_cell_ACRs = DA_cell_types.values()
    combined_ACR_DA_values = pd.concat(all_gathered_cell_ACRs)
    generate_list_ctr_ACRs = combine_ct_ACRs(all_gathered_cell_ACRs)
    
    ## Save at end
    ACR_cell_type_output_name = args.o + ".ACR_celltypes.key.txt"
    generate_list_ctr_ACRs.to_csv(ACR_cell_type_output_name, sep="\t", encoding='utf-8')
    

    celltype_acr_bed_file_output = args.o + ".celltype_limited.bed"
    broadlyaccessible_acr_bed_file_output = args.o + ".broadly_accessible.bed"
    all_acr_output_file_output = args.o + ".all_acrs.bed"

    
    cell_type_bed_file = []
    broadly_accessible_bed_file = []
    all_bed_files = []
        
    ## Moving to dict is like 1000X faster than iloc look ups
    acr_cell_type_dict = dict(zip(generate_list_ctr_ACRs.gene_name, generate_list_ctr_ACRs.cell_type))
    bed_file = read_bed_file(args.bed)
    for line in tqdm(bed_file):
        if line[3] in acr_cell_type_dict:
            ## Assign the cell types in that ACR
            grab_acr_cell_types = acr_cell_type_dict[line[3]]
            final_string = ",".join(grab_acr_cell_types)
            updated_ACR_name = line[3] + ":" + final_string

            ## Generate the new line and add
            new_line = line
            new_line[3] = updated_ACR_name
            cell_type_bed_file.append(new_line)
            all_bed_files.append(new_line)
        elif line[3] not in acr_cell_type_dict:
            broadly_accessible_bed_file.append(line)
            all_bed_files.append(line)

    
    write_output(broadly_accessible_bed_file, broadlyaccessible_acr_bed_file_output)
    write_output(cell_type_bed_file, celltype_acr_bed_file_output)
    write_output(all_bed_files, all_acr_output_file_output)
