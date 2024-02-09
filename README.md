# Scripts Used for the Analysis in the Manuscript Investigating the cis-Regulatory Basis of C3 and C4 Photosynthesis in Grasses at Single-Cell Resolution

Below you will find a detailed list of scripts used in the above manuscript. 
Scripts are docuemnted and ordered by the various languages they've been written in. Below are some specfic entries that may be of interest to some readers.


### QC of sciATAC Libraries and Clustering
All QC and Clustering of sciATAC-seq libraries was done using the software packages [Socrates](https://github.com/plantformatics/Socrates). These steps are outlined on a per-species basis and can be found in 
two nested directories, namely: https://github.com/Jome0169/Mendieta.C4_manuscript/tree/main/ipynbs/QC_libraries and here: https://github.com/Jome0169/Mendieta.C4_manuscript/tree/main/ipynbs/cluster_libraries. 


### Calling Peaks on Annotated Datasets
The peakcalling done in this manuscript is implemented thoroughly in the following snakemake script here: https://github.com/Jome0169/Mendieta.C4_manuscript/tree/main/snakemake/peak_calling. With the important peak calling script being located here:https://github.com/Jome0169/Mendieta.C4_manuscript/blob/main/python_scripts/scATAC/call_scACRs.py. In brief this script takes in a meta file, a Tn5 file of integrations, genome size, and an FDR value in order to identify cell-type-specific ACRs genome wide. An example to run this code is as follows:     

```python
python scripts/call_scACRs.py -bed {input.bed_file} \
    -meta {input.meta_file} \
    -col {params.col_val} \
    -fai {input.fai_file} -bw yes \
    -gsize {params.gsize} \
    -base {params.base_name} -outdir {params.output_dir} -cores \
    {threads} -rep {params.rep_col} -fdr {params.fdr_val}
```
The output from this run will be a set of peaks, as well as a set of normalized bedgraphs, normalized both by the total library per million scaling factor as well as CPM normalized. These bedgraphs can quickly be loaded into a genome browswer and visualzied as such below:
![Screen Shot 2024-02-09 at 11 20 49 AM](https://github.com/Jome0169/Mendieta.C4_manuscript/assets/8882330/9a1e059b-f5c8-4616-b66f-cc3e4cdafc3b)

Output directories will take the following format, with the final list of peaks being named in this example `zm.peaks.500bp_peaks.no_exons.bed`. For each cell-type you'll want to use the bedgraph file with the extension `.normalized.bdg` 

![Screen Shot 2024-02-09 at 11 22 59 AM](https://github.com/Jome0169/Mendieta.C4_manuscript/assets/8882330/f43aff44-a994-4810-9978-85cdcdd0b5f1) 


