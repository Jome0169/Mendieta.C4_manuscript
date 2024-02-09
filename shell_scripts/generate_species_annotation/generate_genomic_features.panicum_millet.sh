#!/bin/bash

# Usage: ./generate_features.sh <input_file> <output_prefix>

input_file="$1"
chrom_size="$2"
output_prefix="$3"

# Generate introns
gt gff3 -retainids -addintrons "$input_file" > "${output_prefix}.gene_exons_introns.gff3"

# Generate exons
bioawk -c gff '$feature == "exon" {print $seqname,$start,$end,$feature,".",$strand}' OFS="\t" "${output_prefix}.gene_exons_introns.gff3" > "${output_prefix}.all_exons.bed"

# Generate introns
bioawk -c gff '$feature == "intron" {print $seqname,$start,$end,$feature,".",$strand}' OFS="\t" "${output_prefix}.gene_exons_introns.gff3" > "${output_prefix}.introns.bed"

# Generate TSS based on "gene" feature, 500 BP upstream of gene
bioawk -c gff '$3=="gene" {if($7=="+") {if ($4-500 >= 1) print $1,$4-500,$4,$feature,".",$7; else print $1,1,$5,$feature,".",$7} else {if ($5+500 >= 1) print $1,$5,$5+500,$feature,".",$7; else print $1,$4,1,$feature,".",$7}}' OFS="\t" "${output_prefix}.gene_exons_introns.gff3" > "${output_prefix}.TSS.bed"

bedtools subtract -a "${output_prefix}.all_exons.bed" -b "${output_prefix}.TSS.bed" > "${output_prefix}.exons.bed"

# Take 2000 BP upstream of TSS to generate promoter
bedtools flank -l 2000 -r 0 -s -i "${output_prefix}.TSS.bed" -g ${chrom_size} | bedtools subtract -a - -b "${output_prefix}.all_exons.bed" "${output_prefix}.TSS.bed" | awk '{print $1,$2,$3,"promoter",".",$6}' OFS="\t" > "${output_prefix}.promoters.bed"

# Generate chromosome length sequence
awk '{print $1 "\t" 1 "\t" $2}' ${chrom_size} > "${output_prefix}.total_chrom.bed"

# Generate intergenic space
bedtools subtract -a "${output_prefix}.total_chrom.bed" -b "${output_prefix}.TSS.bed" "${output_prefix}.promoters.bed" "${output_prefix}.introns.bed" "${output_prefix}.exons.bed" > "${output_prefix}.intergenic_space.bed"

