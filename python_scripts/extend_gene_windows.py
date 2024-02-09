import sys
import pybedtools
import argparse
import os
#from pybedtools import *



def alter_bed_feature(feature, flag_val):
    """
    """


    if flag_val == "up":
        if int(feature[-1]) == -1:
            pass
        else:
            feature.start = int(feature[-1])
        return(feature)

    elif flag_val == "down":
        if int(feature[-1]) == -1:
            pass
        else:
            feature.stop = int(feature[-1])
        return(feature)
    else:
        print("Error: flag_val must be either up or down")
        sys.exit(1)


def update_closest_calls(bed_closest_file, query_genes_index, all_genes_index, up_down):
    """
    """


    group_list = list(range(1, query_genes_index, 1))

    if up_down == "up":
        merged_vals = bed_closest_file.groupby(g=group_list, 
        c=[(query_genes_index + 3 )], o=["min"]).saveas()

        merged_vals = merged_vals.each(alter_bed_feature, flag_val="up").saveas()

        return(merged_vals)

    elif up_down == "down":
        merged_vals =  bed_closest_file.groupby(g=group_list, 
        c=[(query_genes_index + 2 )], o=["max"]).saveas()

        merged_vals = merged_vals.each(alter_bed_feature, flag_val="down").saveas()

        return(merged_vals)


def reformat(x):
    fixed_lines = []
    for line in open(x.fn):
        fields = line.strip().split('\t')
        fixed_line = [fields[0], fields[3], fields[4], fields[1], fields[2]]
        fixed_lines.append(fixed_line)
    return fixed_lines


def combine_bed_final(feature, index_n):
    """
    """

    feature.stop = int(feature[index_n + 2])
    return(feature)

def merge_up_down(upstream_merged, downstream_merged):
    """
    """

    starting_filed_count = upstream_merged.field_count()
    intersecting_values = upstream_merged.intersect(downstream_merged, wa=True, wb=True).each(combine_bed_final, starting_filed_count)
    return(intersecting_values)



def extend_to_n_genes_away(query_genes, all_genes, query_genes_index, all_genes_index, n=2):
    """
    Extend a given gene to n genes away in both directions.
    """


    # Using closest with t="first" to get the upstream gene.
    upstream = query_genes.closest(all_genes, D="ref", id=True, k=n, io = True).filter(lambda x: x[-1] != '0').saveas()
    upstream_merged = update_closest_calls(upstream, query_genes_index, all_genes_index, "up").saveas()

    # Using closest with t="last" to get the downstream gene.
    downstream = query_genes.closest(all_genes, D="ref", iu=True, k=n, io = True).filter(lambda x: x[-1] != '0').saveas()
    downstream_merged = update_closest_calls(downstream, query_genes_index, all_genes_index, "down").saveas()


    #x = upstream_merged.intersect(downstream_merged, wa=True, wb=True).saveas("intersect.bed")
    # generated_final_string = merge_up_down(upstream_merged, downstream_merged).saveas("intersect.bed")
    # combined_values = upstream_merged.cat(downstream_merged, postmerge=False).saveas('extended.bed')
    # final_combined = combined_values.groupby(g=[1,4,5], c = [2,3], o=["min", "max"])
    # generated_final_string = final_combined.saveas("final.bed")


    return(upstream_merged, downstream_merged)

def write_nested_list_to_file(data, filename=None):
    # Convert the nested list to a tab-separated string
    content = '\n'.join(['\t'.join(map(str, row)) for row in data])
    
    # If no filename is provided, print the content
    if not filename:
        print(content)
        return

    # Remove the file if it already exists
    if os.path.exists(filename):
        os.remove(filename)
    
    # Write the content to the file
    with open(filename, 'w') as file:
        file.write(content)

def get_parser():
    parser = argparse.ArgumentParser(description='Finds peaks shared between \
        replicate peak calls, as well as unqiue peaks to each replicate and \
        outputs said peaks. ')
    parser.add_argument('-bed','--bed_file', help='Bed File to Mimic',\
        required=True, dest='bed'),
    #parser.add_argument('-TFs','--TF_file', help='TF file in bed',\
    #    required=True, dest='TF_b'),
    parser.add_argument('-all_genes','--all_genes', help='Bed file of all genes to make bin with ',\
        required=True, dest='all'),
    parser.add_argument('-genome','--genome_file', help='Genome File to use', \
        required=False, dest='gn'),
    parser.add_argument('-genome_index','--genome_index', help='Genome File to use', \
        required=False, dest='gni'),
    parser.add_argument('-n','--number_genes', help='Number of regions to extend up to - defaults is 2', \
        required=False, dest='n'),

    parser.add_argument('-o','--output_name', help='output', \
        required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()

    # Load the genes in BED format
    chosen_regions_bed = pybedtools.BedTool(args.bed).sort().saveas()
    grab_starting_index = chosen_regions_bed.field_count()
    all_gene_in_genomes = pybedtools.BedTool(args.all).sort().saveas()
    all_gene_in_genomes_index = all_gene_in_genomes.field_count()

    
    if args.n == None:
        upstream_regions, downstream_regions = extend_to_n_genes_away(chosen_regions_bed, 
        all_gene_in_genomes, 
        grab_starting_index, 
        all_gene_in_genomes_index, 
        n=2)

    elif args.n != None and int(args.n) >= 1:
        upstream_regions, downstream_regions = extend_to_n_genes_away(chosen_regions_bed, 
        all_gene_in_genomes, 
        grab_starting_index, 
        all_gene_in_genomes_index, 
        n=int(args.n))

    named_region_dictionary = {}
    for i in chosen_regions_bed:
        name = i[3]
        if name not in named_region_dictionary:
            named_region_dictionary[name] = [[],[]]
        else:
            pass

    for region in upstream_regions:
        name = region[3]
        named_region_dictionary[name][0].append(list(region))
    for region in downstream_regions:
        name = region[3]
        named_region_dictionary[name][1].append(list(region))

    updated_bed_regions = []
    for gene_name, regions_IDd in named_region_dictionary.items():
        take_chrom = regions_IDd[0][0][0]
        take_geneId = regions_IDd[0][0][3]
        take_name = regions_IDd[0][0][4]
        take_start = min(int(upstrea_vals[1]) for upstrea_vals in regions_IDd[0])
        take_stop = max(int(upstrea_vals[2]) for upstrea_vals in regions_IDd[1])
        final_bed_line = [take_chrom, str(take_start), str(take_stop), take_geneId, take_name]
        updated_bed_regions.append(final_bed_line)
    
    write_nested_list_to_file(updated_bed_regions, args.o)
    









