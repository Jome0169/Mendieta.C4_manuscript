import sys
import os
from sys import argv

dictionary = {}
with open(argv[1], 'r') as f:
    for line in f:
        cleaned_line = line.strip().split()
        if cleaned_line[4] not in dictionary:
            dictionary[cleaned_line[4]] = [cleaned_line]
        elif cleaned_line[4] in dictionary:
            dictionary[cleaned_line[4]].append(cleaned_line)


for key,val in dictionary.items():
    if len(val) == 1:
        only_list = val[0]
        print('\t'.join(only_list))
    elif len(val) > 1:
        counter = 1
        for item in val:
            new_item = item[4] + "_v4v5count_" + str((counter))
            item[4] = new_item
            print('\t'.join(item))
            counter += 1
