#!/usr/bin/env python
from Bio import SeqIO
from argparse import ArgumentParser


def main():
    print("WORKING")
    args = get_args()
    last_gene = None
    transcripts = []

    for record in SeqIO.parse(args.input, "fasta"):
        gene = record.id.split("_")[0]

        if last_gene is None:
            last_gene = gene

        if last_gene == gene:
            transcripts.append(record)
        else:
            print_longest_transcripts(transcripts, args.number)
            last_gene = gene
            transcripts = [record]

    print_longest_transcripts(transcripts, args.number)


def print_longest_transcripts(transcripts, number):
    longest = sorted(transcripts, key=lambda x: len(x.seq), reverse=True)

    for record in longest[:number]:
        print(">"+record.id, record.seq, sep="\n")


def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="input fasta", required=True)
    parser.add_argument("-n", "--number", help="number of longest transcript per gene", type=int, required=True)
    return parser.parse_args()
