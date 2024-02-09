"""
Purpose of this script is to take move the markers from one species to another
using orthologous sequences. To do so I'm using orthofinder, as well as BlastP.
Additionally this script also takes in the previous BED file and the Zea mays
data. All of this infromation allows me to go from a series of markers in one
species (Mays) and output them to a different species.

Now what we want to do is for every marker gene on the left side find the best
possible match for it in sorghum. This is quite simple. I'm going to load the
list into a dictionary and then filter on the pident, and the escore. From
glancing above this generally gives the best list possible
"""
import argparse
import sys
import os
import copy

def remove_file(file_name):
    """TODO: Docstring for remove_file.
    :returns: TODO

    """

    try:
        os.remove(file_name)
    except OSError:
        pass


def replace_string_sub_list(potential_list, replace_len):
    """TODO: Docstring for replace_string_sub_list.
    :returns: TODO

    """
    if int(replace_len) != 0:
        split_list_on_commas = potential_list.split(",")
        final = [i[:-int(replace_len)] for i in split_list_on_commas]
        return(final)
    elif int(replace_len) == 0:
        split_list_on_commas = potential_list.split(",")
        final = split_list_on_commas
        return(final)



def clean_ortho_names(bed_file, ortho_file):
    """TODO: Docstring for function.
    :returns: TODO

    """


def read_orthofinder_input_db_to_query(ortho_file_name, query_subtstring, db_substring):
    """TODO: Docstring for read_orthofinder_input.
    :returns: TODO

    """
    db_to_query = {}

    with open(ortho_file_name, 'r') as f:
        for line in f:
            cleaned_line = line.strip().replace(" ", "").split('\t')
            if "Orthogroup" in cleaned_line[0]:
                pass
            else:
                db_list = replace_string_sub_list(cleaned_line[2], db_substring)
                query_list = replace_string_sub_list(cleaned_line[1], query_subtstring)
                for db_name in db_list:
                    if db_name not in db_to_query:
                        db_to_query[db_name] = query_list
                    elif db_name not in db_to_query:
                        print("Something weird - gene likely has multiple \
                                Orthologues??? Likely Wrong. Gene ID: %s" % db_name)
                        sys.exit(-1)
    return(db_to_query)


def read_orthofinder_input_query_to_db(ortho_file_name, query_subtstring, db_substring):
    """TODO: Docstring for read_orthofinder_input.
    :returns: TODO

    """
    query_to_db = {}
    with open(ortho_file_name, 'r') as f:
        for line in f:
            cleaned_line = line.strip().replace(" ", "").split('\t')
            if "Orthogroup" in cleaned_line[0]:
                pass
            else:

                db_list = replace_string_sub_list(cleaned_line[2], db_substring)
                query_list = replace_string_sub_list(cleaned_line[1], query_subtstring)

                for query_name in query_list:
                    if query_name not in query_to_db:
                        query_to_db[query_name] = db_list
                    elif query_name not in query_to_db:
                        print("Something weird - gene likely has multiple \
                                Orthologues??? Likely Wrong. Gene ID: %s" % query_name)
                        sys.exit(-1)
    return(query_to_db)



def read_blast_file(input_file):
    gene_hit_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split('\t')
            if cleaned_line[0] not in gene_hit_dict:
                gene_hit_dict[cleaned_line[0]] = [cleaned_line]
            elif cleaned_line[0] in gene_hit_dict:
                gene_hit_dict[cleaned_line[0]].append(cleaned_line)
            else:
                pass
    return(gene_hit_dict)


def count_occurence_query_genes(input_file):
    """TODO: Docstring for count_occurence.

    :arg1: TODO
    :returns: TODO

    """
    gene_hit_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split('\t')
            if cleaned_line[1] not in gene_hit_dict:
                gene_hit_dict[cleaned_line[1]] = [cleaned_line]
            elif cleaned_line[1] in gene_hit_dict:
                gene_hit_dict[cleaned_line[1]].append(cleaned_line)
            else:
                pass
    return(gene_hit_dict)


def filter_blast_dict(blast_dict):
    finalized_dict = {}
    number_query_occurences = {}
    for gene, blast_nested_list in blast_dict.items():
        #https://blog.finxter.com/how-to-find-the-max-of-list-of-lists-in-python/
        take_top_hit = max(blast_nested_list, key = lambda x: float(x[-1]))
        finalized_dict[gene] = take_top_hit

        if take_top_hit[1] not in number_query_occurences:
            number_query_occurences[take_top_hit[1]] = 1
        elif take_top_hit[1] in number_query_occurences:
            number_query_occurences[take_top_hit[1]] += 1

    return (finalized_dict, number_query_occurences)


def read_marker_bed(bed_file):
    """TODO: Docstring for read_marker_bed.

    Input Example:
    chr10   1557308 1561218 Zm00001eb405250 brk3    subsidiary_mother_cell,developing_pavement_cell axillaryBud;leaf
    chr10   10954954        10957199        Zm00001eb408350 ZmAAP6_3        xylem_parenchyma        axillaryBud;leaf;crownRoot;root;tassel;ear
    chr10   28847706        28853887        Zm00001eb410930 omt1    mesophyll       axillaryBud;leaf
    chr10   117179327       117180363       Zm00001eb422000 acl1    bulliform_cell  axillaryBud;leaf
    chr10   117376805       117380718       Zm00001eb422090 ZmSMXL3 protophloem_sieve_element,metaphloem_sieve_element,phloem_sieve_element_precursors,procambial_meristem
    chr10   135702002       135705078       Zm00001eb427470 zyb14   leaf_primordia,leaf_tips        axillaryBud;leaf
    chr10   139559890       139566977       Zm00001eb428740 ocl2    peripheral_zone_SAM,upper_floral_meristem,lower_floral_meristem axillaryBud;leaf;tassel;ear
    chr10   141187279       141196584       Zm00001eb429330 ZmGL3   trichome,atrichoblast   axillaryBud;leaf;crownRoot;root
    chr10   146716167       146717535       Zm00001eb431780 epf1    stomatal_precursor,guard_mother_cell,guard_cell axillaryBud;leaf
    chr10   147566728       147568465       Zm00001eb432140 wox4    procambial_meristem     axillaryBud;leaf;crownRoot;root;tassel;ear
    chr10   149373989       149380403       Zm00001eb433020 ZmMP_1  procambial_meristem     axillaryBud;leaf;crownRoot;root;tassel;ear

    """
    marker_name_dict = {}
    with(open(bed_file, 'r')) as f:
        for line in f:
            cleaned_line = line.strip().split()
            #remove the coordinnates - generate dictionary 
            if cleaned_line[0] == "chr":
                pass
            else:
                generate_updated_list = cleaned_line[3:]
                if generate_updated_list[0] not in marker_name_dict:
                   marker_name_dict[generate_updated_list[0]]  = generate_updated_list 
                elif generate_updated_list[0] in marker_name_dict:
                    print("POTENTIAL ERROR - Gene Name %s in marker list twice???" % generate_updated_list[0])
    return(marker_name_dict)


def read_species_bed(bed_file_name):
    """Takes in the bed file from the species with which markers are
    transferring to. It trims this down to the first four coumns, and then
    uses the BLASTP results to modify the bed file. Final output is a list of
    markers in bed format for downstream processing.
    :Input Example:
        Chr01   22084   23338   Sorbiv5.1_pg612.valid.m1.g      .       +       JGI     gene      .       ID=Sorbiv5.1_pg612.valid.m1.g;Name=Sorbiv5.1_pg612.valid.m1.g
        Chr01   31589   35375   Sobic.001G000200        .       -       JGI     gene    .ID=Sobic.001G000200;Name=Sobic.001G000200
        Chr01   42929   63076   Sorbiv5.1_pg1680.valid.m10.g    .       -       JGI     gene      .       ID=Sorbiv5.1_pg1680.valid.m10.g;Name=Sorbiv5.1_pg1680.valid.m10.g
        Chr01   43889   44492   Sorbiv5.1_pg4994.m1.g   .       -       JGI     gene    .ID=Sorbiv5.1_pg4994.m1.g;Name=Sorbiv5.1_pg4994.m1.g
        Chr01   73507   74211   Sobic.001G000501        .       +       JGI     gene    .ID=Sobic.001G000501;Name=Sobic.001G000501
        Chr01   80923   83931   Sobic.001G000700        .            JGI     gene    .ID=Sobic.001G000700;Name=Sobic.001G000700
        Chr01   83447   90229   Sorbiv5.1_pg4489.valid.m32.g    .       -       JGI     gene      .       ID=Sorbiv5.1_pg4489.valid.m32.g;Name=Sorbiv5.1_pg4489.valid.m32.g
        Chr01   99766   102156  Sobic.001G000900        .       +       JGI     gene    .ID=Sobic.001G000900;Name=Sobic.001G000900+

    :returns:
       Sorbiv5.1_pg612. : [Chr01,22084,23338,Sorbiv5.1_pg612.]
       Sobic.001G000200 : [Chr01,31589,35375,Sobic.001G000200]
       Sorbiv5.1_pg1680 : [Chr01,42929,63076,Sorbiv5.1_pg1680]
       Sorbiv5.1_pg4994 : [Chr01,43889,44492,Sorbiv5.1_pg4994]
       Sobic.001G000501 : [Chr01,73507,74211,Sobic.001G000501]
       Sobic.001G000700 : [Chr01,80923,83931,Sobic.001G000700]
       Sorbiv5.1_pg4489 : [Chr01,83447,90229,Sorbiv5.1_pg4489]
       Sobic.001G000900 : [Chr01,99766,102156 ,Sobic.001G000900]

    """
    gene_bed_dict = {}
    with open(bed_file_name, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split()
            take_first_four = cleaned_line[0:4]
            if take_first_four[3] not in gene_bed_dict:
                gene_bed_dict[take_first_four[3]] = take_first_four
            elif take_first_four[3] in gene_bed_dict:
                gene_bed_dict[take_first_four[3]].append(take_first_four)
    return(gene_bed_dict)


def generate_final_list(blast_dict, marker_dict, species_dict, string_to_add):
    """TODO: Docstring for generate_final_list.
    :returns: TODO

    """
    
    key_list = []
    final_line = []
    for key, val in blast_dict.items():

        replace_marker_name = key.replace("_P001", "")
        pulled_og_marker_info = marker_dict[replace_marker_name]


        #This will likely have to change... Irritating as proteins are 
        #encoded uniqly across species... Also this solution is gross

        fake_protein_range = range(1,100)
        check_isoform_list = [".m" + str(i) for i in fake_protein_range]
        final_tuple = tuple(check_isoform_list)
        print(val[1])
        if val[1].endswith(final_tuple):
            pull_species_bed_line = species_dict[val[1]+".g"]
            final_correct_key = val[1] +'.g'
        elif val[1].startswith("ta") or "promoted" in val[1]: 
            pull_species_bed_line = species_dict[val[1] +'.g']
            final_correct_key = val[1] +'.g'
        else:
            generate_correct_key_for_query  = val[1].split('.')
            final_correct_key = '.'.join(generate_correct_key_for_query[:-1])
            pull_species_bed_line = species_dict[final_correct_key]
        
        # Generate a final marker name so we know what species it originated
        # from 
        generate_updated_marker_line = string_to_add + "_" + pulled_og_marker_info[1]
        additional_marker_info = pulled_og_marker_info[2:]

        query_species_final_line = pull_species_bed_line + [generate_updated_marker_line] + additional_marker_info

        final_line.append(query_species_final_line)


        ## Generate Key File information
        og_species_marker_info = pulled_og_marker_info[0:3]
        marker_key_list = [final_correct_key] + og_species_marker_info 
        key_list.append(marker_key_list)

    return(key_list, final_line)


def merge_duplicates(count_query_dict, final_bed_file, final_key_file):
    """TODO: Docstring for merge_duplicates.
    :returns: TODO

    """


    merge_groups = {}
    merge_keys = {}
    merge_finalized_list = []
    merge_finalized_key_file = []

    all_markers_merged_single = []
    all_key_file_merged_single = []

    not_passing_gene_names = []
    for gene_name_query, count in count_query_dict.items():
        for line in final_bed_file:
            
            if line[3].endswith(".g"):
                gene_name_bed = line[3].replace(".g", "")
            else:
                gene_name_bed = line[3]
        
            if gene_name_bed in gene_name_query and gene_name_query not in merge_groups and count >1 :
                merge_groups[gene_name_query] = [line]
                not_passing_gene_names.append(gene_name_bed)
            elif gene_name_bed in gene_name_query and gene_name_query in merge_groups and count >1:
                merge_groups[gene_name_query].append(line)
                not_passing_gene_names.append(gene_name_bed)
            #elif count == 1:
            #    pass
            #    all_markers_merged_single.append(line)

        for key_line in final_key_file:
            if key_line[0].endswith(".g"):
                gene_name_bed = key_line[0].replace(".g", "")
            else:
                gene_name_bed = key_line[0]
            if gene_name_bed in gene_name_query and count >1 :
                if gene_name_bed not in merge_keys:
                    merge_keys[key_line[0]] = [key_line]
                    not_passing_gene_names.append(gene_name_bed)
                elif gene_name_bed in merge_keys and count >1:
                    merge_keys[key_line[0]].append(key_line)
                    not_passing_gene_names.append(gene_name_bed)
            #else:
            #    all_key_file_merged_single.append(key_line)
    
    for line in final_bed_file:
        if line[3].endswith(".g"):
            gene_name_bed = line[3].replace(".g", "")
        else:
            gene_name_bed = line[3]
        if gene_name_bed not in not_passing_gene_names:
            all_markers_merged_single.append(line)

    for line in final_key_file:
        if line[0].endswith(".g"):
            gene_name_bed = line[0].replace(".g", "")
        else:
            gene_name_bed = line[0]
        if gene_name_bed not in not_passing_gene_names:
            all_key_file_merged_single.append(line)


    


    non_merged_bed_lines = [list(y) for y in set([tuple(x) for x in all_markers_merged_single])]
    non_merged_key_files = [list(y) for y in set([tuple(x) for x in all_key_file_merged_single ])]

    for key, val in merge_groups.items():

        grab_same_segment = val[0][0:4]
        merged_marker_name = [i[4] for i in val]
        joined_marker_names = '__'.join(merged_marker_name)
        cell_type = [i[5].split(",") for i in val]
        flat_cell_information = ",".join(set([item for sublist in cell_type for item in sublist]))

        tissue_type = [i[6].split(";") for i in val]
        flat_tissue_type = ";".join(set([item for sublist in tissue_type for item in sublist]))

        final_bed_line = grab_same_segment + [joined_marker_names] + [flat_cell_information] + [flat_tissue_type]
        merge_finalized_list.append(final_bed_line)


    for key, val in merge_keys.items():

        combined_gene_names = "__".join(z[1] for z in val)
        combined_marker_names = "__".join(z[2] for z in val)
        combined_cell_information = [z[3].split(',') for z in val]
        flat_cell_information = ",".join(set([item for sublist in combined_cell_information for item in sublist]))

        final_line = [key, combined_gene_names, combined_marker_names, flat_cell_information]
        merge_finalized_key_file.append(final_line)

    
    final_bed_lines = non_merged_bed_lines + merge_finalized_list
    final_key_file = non_merged_key_files + merge_finalized_key_file



    return(final_bed_lines, final_key_file)



def get_parser():
    parser = argparse.ArgumentParser(description='generated a sparse matrix for \
            later analysis')
    
    parser.add_argument('-marker_bed', '--marker_bed',  
        required=True, dest ="marker_bd"),
    parser.add_argument('-len_sub_query', '--len_sub_query',  
        required=True, dest ="len_sub_query"),
    parser.add_argument('-len_sub_df', '--len_sub_df',  
        required=True, dest ="len_sub_df"),
    parser.add_argument('-ortho_query', '--ortho_query',  
        required=True, dest ="ortho_query"),
    parser.add_argument('-query_bed', '--query_bed',  
        required=True, dest ="query_bed"),
    parser.add_argument('-o', '--output_base',
           required=True, dest ="output")
    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()


    #Read in orthofinder results:
    db_to_query_ortho_dict = read_orthofinder_input_db_to_query(args.ortho_query, args.len_sub_query, args.len_sub_df)
    #query_to_db_ortho_dict = read_orthofinder_input_query_to_db(args.ortho_query, args.len_sub_query, args.len_sub_df)


    #ref_db_to_query_ortho_dict = read_orthofinder_input_db_to_query(args.ortho_ref, args.len_sub_df, args.len_sub_query )
    #ref_query_to_db_ortho_dict = read_orthofinder_input_query_to_db(args.ortho_ref, args.len_sub_df, args.len_sub_query)
   
    #Parse  marker dict
    marker_dict = read_marker_bed(args.marker_bd)
    species_bed_file = read_species_bed(args.query_bed)


    final_marker_list_bed = []
    final_marker_key_list = []
    final_marker_dict = {}
    no_markers_present = []
    counter_of_query = {}

    for key, value in marker_dict.items():
        protein_name = key
        if protein_name in db_to_query_ortho_dict: 


            grabbed_ortho = db_to_query_ortho_dict[protein_name]
            
            if len(grabbed_ortho) == 1:

                grabbed_ortho = grabbed_ortho[0]
                #Match the protein sequence from orthofinder with bed file by
                #matching substrings in the key file
                grab_ortho_bed_name = [x for x in species_bed_file.keys() if grabbed_ortho in x]
                                                                                                 
                if len(grab_ortho_bed_name) > 0:
                    #use the above substring to pull the bed file
                    grabbed_bed_file = species_bed_file[grab_ortho_bed_name[0]]
                    #Grab the reference marker information
                    grabbed_reference_marker = marker_dict[protein_name]
                    
                    #The 1+ Regions is because we don't want the maize gene name
                    generate_combined_list = grabbed_bed_file + grabbed_reference_marker[1:]
                
                
                    #This way if two db markers hit to the same gene, we can fix it
                    #in the final analysis of this dictionary
                    if grab_ortho_bed_name[0] not in final_marker_dict:
                        final_marker_dict[grab_ortho_bed_name[0]] = [generate_combined_list]
                    elif grab_ortho_bed_name[0] in final_marker_dict:
                        final_marker_dict[grab_ortho_bed_name[0]].append(generate_combined_list)
                                                                                                     
                    #Generate a sub list for the key file
                    key_list = grab_ortho_bed_name + grabbed_reference_marker
                    final_marker_key_list.append(key_list)
                else:
                    pass
            elif len(grabbed_ortho) > 1 and len(grabbed_ortho) != 0:
                grab_additonal_ortho = db_to_query_ortho_dict[protein_name]
                take_number_orthos = range(1,len(grabbed_ortho) + 1)

                for extra_marker, number in zip(grabbed_ortho, take_number_orthos):
                    #Match the protein sequence from orthofinder with bed file by
                    #matching substrings in the key file
                    grab_ortho_bed_name = [x for x in species_bed_file.keys() if extra_marker in x]

                    #use the above substring to pull the bed file
                    try:
                        grabbed_bed_file = species_bed_file[grab_ortho_bed_name[0]]
                        #Grab the reference marker information
                        grabbed_reference_marker = marker_dict[protein_name]

                        number_reference_list_marker = grabbed_reference_marker[1] + "_ortho_count_" + str(number)
                        updated_list_marker = copy.deepcopy(grabbed_reference_marker)
                        updated_list_marker[1] = number_reference_list_marker
                        #The 1+ Regions is because we don't want the maize gene name
                        generate_combined_list = grabbed_bed_file + updated_list_marker[1:]
                        #final_marker_list_bed.append(generate_combined_list)

                        if grab_ortho_bed_name[0] not in final_marker_dict:
                            final_marker_dict[grab_ortho_bed_name[0]] = [generate_combined_list]
                        elif grab_ortho_bed_name[0] in final_marker_dict:
                            final_marker_dict[grab_ortho_bed_name[0]].append(generate_combined_list)
                        
                        #Generate a sub list for the key file
                        key_list = grab_ortho_bed_name + grabbed_reference_marker
                        final_marker_key_list.append(key_list)
                    except:
                        pass
        elif protein_name not in db_to_query_ortho_dict:

            grabbed_reference_marker = marker_dict[protein_name]
            no_markers_present.append(grabbed_reference_marker)

           

    for query_marker, list_vals in final_marker_dict.items():
        if len(list_vals) == 1:
            final_marker_list_bed.append(list_vals[0])
        elif len(list_vals) > 1:
            base_list = list_vals[0][0:4]
            all_markers = []
            all_marker_cell_types = []
            all_marker_tissues = []

            for sub_list in list_vals:
                print(sub_list)
                all_markers.append(sub_list[4])
                all_marker_cell_types.append(sub_list[5])
                all_marker_tissues.append(sub_list[6])
            
            joined_marker_names = "__".join(all_markers)
            joined_marker_cell_tupes = ",".join(set(all_marker_cell_types))
            joined_marker_tissues = ";".join(set(all_marker_tissues))
            final_new_marker_bed =  base_list + [joined_marker_names] + [joined_marker_cell_tupes] + [joined_marker_tissues]
            final_marker_list_bed.append(final_new_marker_bed)

    key_file_name = args.output + '.key_file.txt'
    marker_name = args.output + '.markers.bed'
    non_markers_found = args.output + '.markers_not_found.txt'

    remove_file(key_file_name)
    remove_file(marker_name)
    remove_file(non_markers_found)

    with open(key_file_name, 'a+') as f:
        for list_n in final_marker_key_list:
            f.write('\t'.join(list_n))
            f.write('\n')

    with open(marker_name, 'a+') as f:
        for list_n in final_marker_list_bed:
            f.write('\t'.join(list_n))
            f.write('\n')                   

    with open(non_markers_found, 'a+') as f:
        for list_n in no_markers_present:
            f.write('\t'.join(list_n))
            f.write('\n')                   


