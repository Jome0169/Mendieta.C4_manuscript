
import argparse
import sys
import os
import pysam 
import pybedtools

def read_fasta_file(fasta_file):
    """Takes in bed file and reads it in using pybedtools. Will make life
    easier longitudinally.

    :bed_file: TODO
    :returns: TODO

    """
    try:
        os.path.isfile(bed_file)
    except:
        FileNotFoundError
    generate_header_seq_dict = {}
    with open(fasta_file, 'r') as f:
         for line in f:
            if line.startswith(">"):
                 header_line = line.strip()
            else:
                sequence_line = line.strip()
                generate_header_seq_dict[header_line] = sequence_line
    return(generate_header_seq_dict)

def parse_fasta_bed(fasta_dict):
    generated_cell_type_nested_dict = {}

    for key, val in fasta_dict.items():
        generated_split_key = key.split("::")[0]
        take_acr_type = generated_split_key.split(";")[1]
        if "," in take_acr_type:
            pass
        elif "," not in generated_cell_type_nested_dict:
            if take_acr_type not in generated_cell_type_nested_dict:
                generated_cell_type_nested_dict[take_acr_type] = {key: val}
            elif take_acr_type in generated_cell_type_nested_dict:
                generated_cell_type_nested_dict[take_acr_type][key] = val

    return generated_cell_type_nested_dict



def write_output_fastas(base_name, cell_type_dict):
    for cell_type, seqs in cell_type_dict.items():
        output_file_name = f"{base_name}.acrs_{cell_type}.fa"
        with open(output_file_name, 'a') as f:
            for name, seq in seqs.items():
                f.write(name)
                f.write("\n")
                f.write(seq)
                f.write("\n")



def get_parser():
    parser = argparse.ArgumentParser(description='Calculate regions with fold \
            enrichment over X and outputs these regions.')
    parser.add_argument('-fasta', "--fa", help="input file", required=True, \
                        dest='fa')
    parser.add_argument('-o', "--output", help="outptu base file name. If not \
    given will be print ", required=True, dest='o')
    args = vars(parser.parse_args())    
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()


    print("Working on %s" % args.fa)
    fa_seqs = read_fasta_file(args.fa)
    cell_types = parse_fasta_bed(fa_seqs)
    write_output_fastas(args.o, cell_types)


    

