import argparse
import sys
import os


"""
Parses output from find_cns regions. Takes file input like:
#cns_id,qaccn,qchr,qstart,qstop,qstrand,saccn,schr,sstart,sstop,sstrand,eval,postion,qfasta,sfasta,link

Sb007G166800|7|60179223|60179237|+|Zm00001d031911|1|206536393|206536407|+|29.5,Sb007G166800,7,60179223,60179237,+,Zm00001d031911,1,206536393,206536407,+,29.5,intron,TTGAACATTTTTCA,TTGAACATTTTTCA,http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&show_cns=1&dsgid1=28853&chr1=7&x1=60179223&dr1up=15000&dr1down=15000&dsgid2=33766;chr2=1;x2=206536393;dr2up=15000;dr2down=15000;num_seqs=2;hsp_overlap_limit=0;hsp_size_limit=0

and generates bed files for the two species. Here being zea mays, and sorghum.
HArdcoded index pulls are below.
"""


def read_file(input_file):
    lines_parsed = []
    with open(input_file, 'r') as f:
        for line in f:
            cleaned_line = line.strip().split(',')
            lines_parsed.append(cleaned_line)
    return(lines_parsed)

def parse_lines(list_nested_lines):
    species_1_list = []
    species_2_list = []

    for i in list_nested_lines:
        if i[0].startswith('#'):
            pass
        else:
            species_one_values = [i[2], i[3], i[4],i[1],".",i[5]]
            species_two_values = [i[7], i[8], i[9], i[6],".",i[10]]

            species_1_list.append(species_one_values)
            species_2_list.append(species_two_values)
    return(species_1_list, species_2_list)


def write_bed_file(speces_name, nested_list):

    generate_output_file = speces_name + ".CNS.bed"
    with open(generate_output_file, 'a') as f:
        for i in nested_list:
            f.write('\t'.join(i))
            f.write('\n')





def get_parser():
    parser = argparse.ArgumentParser(description='generated a sparse matrix for \
            later analysis')
    parser.add_argument('-i', '--input-file',  
            required=True, dest ="input")
    parser.add_argument('-sp1', '--species1',
           required=True, dest ="sp1")
    parser.add_argument('-sp2', '--species2',
           required=True, dest ="sp2")


    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()
    
    CNS_list = read_file(args.input)
    species_1_CNS, species_2_CNS = parse_lines(CNS_list)

    write_bed_file(args.sp1, species_1_CNS)
    write_bed_file(args.sp2, species_2_CNS)









