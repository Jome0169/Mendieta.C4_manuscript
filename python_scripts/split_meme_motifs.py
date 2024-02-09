"""
split_meme_motifs.py 

Take in a .meme file with the basis looking like:

MEME version 4.4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from /mnt/thumper-e1/home/shhuang/projects/dap/analysis.v4/gem07_rep_memechip03/ABI3VP1_tnt/AT5G18090_col_a/background):
A 0.30000 C 0.20000 G 0.20000 T 0.30000

MOTIF ABI3VP1_tnt.AT5G18090_col_a_m1 AT5G18090

letter-probability matrix: alength= 4 w= 15 nsites= 147 E= 2.1e-117
  0.244898        0.163265        0.401361        0.190476
  0.414966        0.149660        0.238095        0.197279
  0.312925        0.176871        0.129252        0.380952
  0.000000        0.000000        1.000000        0.000000
  1.000000        0.000000        0.000000        0.000000
  0.000000        0.000000        0.000000        1.000000
  0.000000        0.006803        0.993197        0.000000
  1.000000        0.000000        0.000000        0.000000
  1.000000        0.000000        0.000000        0.000000
  0.292517        0.122449        0.360544        0.224490
  0.394558        0.258503        0.176871        0.170068
  0.380952        0.190476        0.204082        0.224490
  0.340136        0.163265        0.278912        0.217687
  0.435374        0.142857        0.210884        0.210884
  0.265306        0.102041        0.326531        0.306122

URL http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php?AGI=AT5G18090


MOTIF ABI3VP1_tnt.AT5G25475_col_a_m1 AT5G25475

letter-probability matrix: alength= 4 w= 6 nsites= 218 E= 2.3e-147
  0.000000        0.899083        0.100917        0.000000
  1.000000        0.000000        0.000000        0.000000
  1.000000        0.000000        0.000000        0.000000
  0.000000        0.000000        1.000000        0.000000
  0.000000        1.000000        0.000000        0.000000
  1.000000        0.000000        0.000000        0.000000

URL http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php?AGI=AT5G25475


And split this into the corresponding motfis with a single motif being in each
file. 
"""

import argparse
import sys
import os

def remove_if_exists(filename):
    """TODO: Docstring for remove_if_exists.
    :returns: TODO
    """
    try:
        os.remove(filename)
    except OSError:
        pass

def read_file(arg_file):
    """
    """

    header_info = []
    motif_info =  []
    
    with open(arg_file, "r") as f:
        meme_by_meme_list = []
        status = 0
        for line in f:
            clean_line = line.strip()
            clean_line_split = line.strip().split()
            if clean_line == None:
                pass

            elif clean_line.startswith("MEME") or \
            clean_line.startswith("ALPHABET") or \
            clean_line.startswith("strands") or \
            clean_line.startswith("Background") or \
            clean_line.startswith("A"):
            

                header_info.append(clean_line)

            elif clean_line.startswith("MOTIF") and status == 0:
                meme_by_meme_list.append([clean_line])
                status = 1

            elif not clean_line.startswith("MOTIF") and clean_line.startswith("letter") and status == 1:
                meme_by_meme_list.append([clean_line])

            elif not clean_line.startswith("MOTIF") and not clean_line.startswith("letter") and status == 1:
                meme_by_meme_list.append(clean_line_split)

            elif clean_line.startswith("MOTIF") and status == 1:
                motif_info.append(meme_by_meme_list)
                meme_by_meme_list = []
                meme_by_meme_list.append(clean_line_split)
            else:
                pass
    return(header_info, motif_info)


def generate_motif_file_name(motif_info_list):
    take_motif_name = motif_info_list[0]    
    no_base_motif = take_motif_name[1:]
    replaced_spaced = '_'.join(no_base_motif).replace(" ", "_") + ".meme"
    return replaced_spaced 


def generte_seperate_file(header_info, motif_info, output_arg):
    for motif_present in motif_info:
        motif_file_name = generate_motif_file_name(motif_present)

        if output_arg != None:
            new_motif_file_name = output_arg + "/" +motif_file_name
        elif output_arg == None:
            new_motif_file_name = motif_file_name

        remove_if_exists(new_motif_file_name)
        with open(new_motif_file_name, 'w+') as w_file:
            for line in header_info:
                if line.startswith("Background "):
                    w_file.write(line)
                    w_file.write('\n')
                else:
                    w_file.write(line)
                    w_file.write('\n')
                    w_file.write('\n')
            for motif_line in motif_present:
                if motif_line != None and len(motif_line) == 1:
                    if motif_line[0].startswith("MOTIF"):
                        w_file.write(motif_line[0])
                        w_file.write('\n')
                        w_file.write('\n')
                    else:
                        w_file.write(motif_line[0])
                        w_file.write('\n')


                elif motif_line != None and len(motif_line) > 1:
                    w_file.write('\t'.join(motif_line))
                    w_file.write('\n')


def get_parser():
    parser = argparse.ArgumentParser(description='generated a sparse matrix for \
            later analysis')
    parser.add_argument('-i', '--input-file',  
            required=True, dest ="input")
    parser.add_argument('-o', '--output_base',
           required=False, dest ="output")
    args = vars(parser.parse_args())
    return parser



if __name__ == "__main__":
    args = get_parser().parse_args()

    head_info, motif_info = read_file(args.input)
    
    generte_seperate_file(head_info, motif_info, args.output)
    
    

