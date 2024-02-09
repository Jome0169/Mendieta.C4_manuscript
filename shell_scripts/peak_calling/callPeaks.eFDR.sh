#!/bin/bash

#PBS -S /bin/bash                       ## shell selection
#PBS -q batch                           ## queue selection
#PBS -N j_callPeaks                     ## job name
#PBS -l nodes=1:ppn=10                  ## ppn=threads per node, Needs to match the software argument integer
#PBS -l walltime=00:12:00:00            ## total time running limit
#PBS -l mem=50gb                       ## memory limit

# set env
#cd $PBS_O_WORKDIR
#source ~/.zshrc

# load modules
#ml MACS2/2.1.1.20160309-foss-2016b-Python-2.7.14
#ml parallel/20160622-foss-2016b

set -euxo pipefail 

# vars
tissue='Athali_root'
threads=10
clusters=Thali_root_merged_final.LST.clusters.fixed.txt
bam=Arabidopsis_tis_Root_all_reps_merged.TN5.no_org.bed
perl -ne 'chomp;if($_=~/^cell/){next;}else{print"$_\n";}' $clusters > $clusters.new
ref=Athal.chrom.sizes
input=/scratch/apm25309/single_cell/ATACseq/v3/bedfiles/tn5_center_50bp/B73input.Tn5.clean.bed

# make directory
if [ ! -d bigwigs ];then
	mkdir bigwigs
else
	rm -rf bigwigs
	mkdir bigwigs
fi

if [ ! -d bedgraphs ]; then
	mkdir bedgraphs
else
	rm -rf bedgraphs
	mkdir bedgraphs
fi

#################################################
# function to merge cells given cluster and BAM #
#################################################
mergeCells(){
	
	# vars
	clust=$1
	bamf=$2
	clustf=$3
	tissue=$4

	# prep cluster data file
	perl -ne 'chomp;if($_=~/^cellID/){next;}else{print"$_\n";}' $clustf > $clustf.$clust
	cut -f1,9 $clustf.$clust > $clustf.$clust.temp
	mv $clustf.$clust.temp $clustf.$clust
	
	# select cells from correct cluster
	awk -F'\t' -v cluster=$clust '$2==cluster' $clustf.$clust | cut -f1 - > cluster.$clust.bc_IDs.txt
	grep -Ff cluster.$clust.bc_IDs.txt $bamf | cut -f1-4 - > cluster_$clust.$tissue.pool.bed

	# make pseudoreps
	numcells=$( wc -l < cluster.$clust.bc_IDs.txt )
	repnum=$( echo "$numcells / 2" | bc -l )
	repnum=$(printf %.0f $repnum)
	echo " - splitting cluster into 2 reps each with $repnum cells ..."
	shuf cluster.$clust.bc_IDs.txt > cluster.$clust.bc_IDs.shuffled.txt
	head -n $repnum cluster.$clust.bc_IDs.shuffled.txt > cluster.$clust.bc_IDs.rep1.txt
	tail -n $repnum cluster.$clust.bc_IDs.shuffled.txt > cluster.$clust.bc_IDs.rep2.txt
	grep -Ff cluster.$clust.bc_IDs.rep1.txt cluster_$clust.$tissue.pool.bed > cluster_$clust.$tissue.rep1.bed	
	grep -Ff cluster.$clust.bc_IDs.rep2.txt cluster_$clust.$tissue.pool.bed > cluster_$clust.$tissue.rep2.bed

	## clean temp files
	rm cluster.$clust.bc_IDs.txt
	rm $clustf.$clust
	rm cluster.$clust.bc_IDs.shuffled.txt
	rm cluster.$clust.bc_IDs.rep1.txt
	rm cluster.$clust.bc_IDs.rep2.txt
	
}
export -f mergeCells


#############################
# iterate over all clusters #
#############################

# function
iterateClusters(){

	# load parameters
	i=$1
	bam=$2
	clusters=$3
	tissue=$4
	ref=$5
	input=$6

	# make directory
	if [ ! -d $PWD/$tissue.$i ]; then
		mkdir $PWD/$tissue.$i
	fi

	# verbose
	echo "merging cells from cluster $i ..."
	mergeCells $i $bam $clusters $tissue

	Done# call peaks
	echo "calling peaks with macs2 cluster $i ..."
	readdepth=$( wc -l < cluster_$i.$tissue.pool.bed)
	echo "total reads = $readdepth"
	macs2 callpeak -t cluster_$i.$tissue.pool.bed \
		-f BED \
        	-g 1.6e9 \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -50 \
		--qvalue 0.05 \
                --outdir $PWD/$tissue.$i \
		--bdg \
                -n cluster.$i.macs2
	
	# un-corrected
	sort -k1,1 -k2,2n $PWD/$tissue.$i/cluster.$i.macs2_treat_pileup.bdg | perl cleanBED.pl $ref $readdepth - > cluster.$i.macs2_treat_pileup.clean.bdg
	bedGraphToBigWig cluster.$i.macs2_treat_pileup.clean.bdg $ref cluster_$i.$tissue.raw.bw
	mv cluster.$i.macs2_treat_pileup.clean.bdg cluster_$i.$tissue.raw.bdg

	# clean
	rm $PWD/$tissue.$i/cluster.$i.macs2_treat_pileup.bdg
	rm $PWD/$tissue.$i/cluster.$i.macs2_control_lambda.bdg


	# call peaks rep1
        echo "calling peaks with macs2 cluster $i ..."
        readdepth=$( wc -l < cluster_$i.$tissue.rep1.bed)
        echo "total reads rep1 = $readdepth"
        macs2 callpeak -t cluster_$i.$tissue.rep1.bed \
                -f BED \
                -g 1.6e9 \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -50 \
                --qvalue 0.05 \
                --outdir $PWD/$tissue.$i \
                -n cluster.$i.macs2.rep1

	# call peaks rep2
        echo "calling peaks with macs2 cluster $i ..."
        readdepth=$( wc -l < cluster_$i.$tissue.rep2.bed)
        echo "total reads rep2 = $readdepth"
        macs2 callpeak -t cluster_$i.$tissue.rep2.bed \
                -f BED \
                -g 1.6e9 \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -50 \
                --qvalue 0.05 \
                --outdir $PWD/$tissue.$i \
                -n cluster.$i.macs2.rep2

	# check overlaps
	bedtools intersect -a $PWD/$tissue.$i/cluster.$i.macs2_peaks.narrowPeak -b $PWD/$tissue.$i/cluster.$i.macs2.rep1_peaks.narrowPeak -u -f 0.1 \
		| bedtools intersect -a - -b $PWD/$tissue.$i/cluster.$i.macs2.rep2_peaks.narrowPeak -u -f 0.1 > $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak


    #Shuffle reproducible peaks
    echo "Generating Shuffled Null Peaks"
    bedtools shuffle -i $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak \
    -excl $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak \
    -g ${ref} \
    | sort -k1,1 -k2,2n - > $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.control

    echo "Counting integration events in Null peaks"
  # count Tn5 sites per permuted region
    bedtools intersect -a $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.control -b cluster_$i.$tissue.pool.bed -c -sorted \
    > $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.control.tn5 

    echo "Counting integration events in Real peaks"
    bedtools intersect -a $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak -b cluster_$i.$tissue.pool.bed -c -sorted \
    > $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.tn5 

    # run eFDR
    echo "FIltering real peaks by FDR value .0001"
    Rscript eFDR_filter_ACRs.R $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.tn5 \
    $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.control.tn5 \
    0.001 \
    $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.final


    echo "Grabbing Summits"
	# get summits for reproducible peaks
	bedtools intersect -a $PWD/$tissue.$i/cluster.$i.macs2_summits.bed \
    -b $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.final -u > $PWD/cluster.$i.reproducible_summits.bed

    echo "Moving Files"
	# move bed files to cluster directory
	mv $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak.final cluster_$i.$tissue.pool.bed cluster_$i.$tissue.rep1.bed cluster_$i.$tissue.rep2.bed $PWD/$tissue.$i


    echo "Done"

}
export -f iterateClusters

# iterate over clusters
parallel -j $threads iterateClusters {1} $bam $clusters.new $tissue $ref $input ::: $( cut -f9 $clusters.new | sort -k1,1n | uniq )

# remove temp1
rm $clusters.new

# create merged set of peaks
./adjustPeaks.sh $tissue




#  # remove blacklist overlaps
#  bedtools intersect -wa -v \
#  -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.narrowPeak \
#  -b ${blacklist} \
#  | sort -k1,1 -k2,2n - > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak
#
#
#  ################################
#  ## Do empirical FDR filtering ##
#  ################################
#
#  # count tn5 sites per ACR
#  bedtools bamtobed -i ${tbam_home_dir}/${INPUT_GENOME}.allreps_${RECIPIENT_GENOME}_no-dups-mapq30.bam \
#    | uniq - \
#    | bedtools intersect -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak -b - -c \
#    > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.tn5
#
#  # permute random regions (excluding ACRs)
#  bedtools shuffle -i ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak \
#    -excl ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.narrowPeak \
#    -g ${contigs_no_genes} \
#    | sort -k1,1 -k2,2n - > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control
#
  # count Tn5 sites per permuted region
#  bedtools bamtobed -i ${tbam_home_dir}/${INPUT_GENOME}.allreps_${RECIPIENT_GENOME}_no-dups-mapq30.bam \
#    | uniq - \
#    | bedtools intersect -a ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control -b - -c \
#    > ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control.tn5
#
#  # run eFDR
#  Rscript eFDR_filter_ACRs.R ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.tn5 \
#    ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control.tn5 \
#    ${eFDR_threshold} \
#    ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak

#  # remove intermediate files
#  rm ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.tn5
#  rm ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control.tn5
#  rm ${out_dir}/${INPUT_GENOME}.${RECIPIENT_GENOME}_peaks.no_blacklist.narrowPeak.control
#
#
#
