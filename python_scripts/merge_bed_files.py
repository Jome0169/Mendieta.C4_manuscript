import argparse
import itertools
import os



def read_files(file_name, known_loci):
    cleaned_lines = []
    with open(file_name, 'r+') as f:
        for line in f:
            clean_line = line.strip().split("\t")
            if clean_line[4] in known_loci:
                pass
            elif clean_line[4] not in known_loci:
                cleaned_lines.append(clean_line)
                known_loci.append(clean_line[4])
    return(cleaned_lines, known_loci)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine bed files.')
    parser.add_argument('-f', '--files', nargs='+', help='<Required> Set flag', required=True)
    args = parser.parse_args()

    k_loci = []
    nested_list = []
    for i in args.files:
        bed_file, k_loci = read_files(i, k_loci)
        nested_list.append(bed_file)

    flattened_list = list(itertools.chain.from_iterable(nested_list))

    for i in flattened_list:
        print("\t".join(i))
