#!/usr/bin/env python

##this script transfers perl script from Alex to python format

#script Updated Mon Jan  4 13:23:43 EST 2021 by Pablo Mendieta
#Script rewrite intially done by Hosung. Currently updating to iter by line and
#write,  as well as take a possibly output file name. This will avoid us
#loading files into memory, and forgo the possability of at segfault. Also
#additional pliability.

import argparse
import re
import sys
import numpy as np
import os


def makeTn5bed(input_sam_fl,output_dir):

    ##transfer the
    ##prepare the flag_dic
    flag_dic = {'0':'read_paired',
                '1':'read_properly',
                '2':'read_umapped',
                '3':'mate_unmapped',
                '4':'read_reverse',
                '5':'mate_reverse',
                '6':'first_pair',
                '7':'second_pair',
                '8':'secondary',
                '9':'fail_qc',
                '10':'duplicate',
                '11':'sup_align'
    }

    store_final_line_list = []

    #with open (input_sam_fl, 'r') as ipt:

    for eachline in input_sam_fl:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        ##extracct barcode information
        bc = ''
        for i in range(9,len(col)):
            if col[i].startswith('BC:Z:'):
                bc = col[i]

        ##store flag list information
        flag_list = [] ##[read_reverse,second_pair,duplicate,sup_align]

        ##transfer the flag to binary information
        bin = np.binary_repr(int(col[1]), width=12)
    

        bin_list = list(bin)
        ##generate flag information
        for i in range(len(bin_list)):
            if bin_list[i] == '1':
                flag_list.append(flag_dic[str(i)])


        ##get start and end pos of read
        chr = col[2]
        pos1 = col[3]
        cigar = col[5]
        netdif = 1

        ##in order to make act list we need to do transfer the CIGAR string to another format:
        ##dic_list is eg. [{'M': 76}, {'I': 15}, {'M': 57}, {'S': 3}]
        list_CIGAR = re.findall('\d+|\D+', cigar)
        n = int(len(list_CIGAR) / 2)
        dic_list = []
        for i in range(0, int(n)):
            dic = {list_CIGAR[(2 * i + 1)]: int(list_CIGAR[(2 * i)])}
            dic_list.append(dic)
        
        cig_dic = {} ####cig_dic stores key and value. key is the 1_M and value is number besides the key in the CIGAR eg {'1_M':50,'2_S':30}
        act_list = [] ##transfer the cigar to [1_M,2_S] or others [1_S,2_M,3_S] stored in the act_list
        ##generate act_list
        item_count = 0
        for eachdic in dic_list:
            item_count += 1
            item_str = str(item_count) + '_' + list(eachdic.keys())[0]
            act_list.append(item_str)
            cig_dic[item_str] = str(eachdic[list(eachdic.keys())[0]])
       
        for eachact in act_list:
            ##eachact is 1_M or 2_S or others
            ##eachact_list = [1,M] or [2,S]
            ##eachact_list[1] is 'M'
            eachact_list = eachact.split('_')

            ##save the time value to the cig dictionary
            ##time indicates the value in the CIGAR eg. 50M. time is '50'
            time = cig_dic[eachact]

            if eachact_list[1] == 'M' or eachact_list[1] == 'S':
                netdif += int(time)

            elif eachact_list[1] == 'D':
                netdif += int(time)
        
        #print("End Net difference")
        #print(netdif) 
        #print("Updated position")
        pos2 = netdif + int(pos1)
        
        ##shift
        if flag_list[0] == 'read_reverse':
            end = pos2 - 4
            start = end - 1
            #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-')
            final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '-'

            if output_dir != None:
                store_final_line_list.append(final_line)
            else:
                print(final_line)

        else:
            start = int(pos1) + 5
            end = start + 1
            #print (chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+')
            final_line = chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(bc) + '\t' + '+'
            if output_dir != None:
                store_final_line_list.append(final_line)
            else:
                print(final_line)


    if output_dir != None:
        with open(output_dir, 'a') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

    elif output_dir == None:
        pass




def remove_output_file(file_name):
    """TODO: Docstring for remove_output_file.
    :returns: TODO

    """
    try:
        os.remove(file_name)
    except OSError:
        pass



def get_parser():
    parser = argparse.ArgumentParser(description='Pull our reads aligning to a\
        region from multiple list of BAM files, puts them into a BAM file\
        for later assembly.')
    parser.add_argument('-sam', "--sam_file", help="Bam file to \
        pull reads from.", required=True, dest='sam', nargs="?", \
        type = argparse.FileType("r"), default=sys.stdin)
    parser.add_argument('-output_file', "--output", help="Output file to write to. \
            If none given output writes to sOutput file to write to. If none \
            given output writes to sout.", required=False, dest='o')

    args = vars(parser.parse_args())
    return parser


if __name__ == "__main__":
    args = get_parser().parse_args()
    
    #Load all Bed files

    #final_lines = makeTn5bed(args.bam)

    if args.o != None:
        remove_output_file(args.o)
        makeTn5bed(args.sam,args.o)
    if args.o == None:
        makeTn5bed(args.sam, None)
