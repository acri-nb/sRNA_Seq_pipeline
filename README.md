# Derfinder_pipeline implementation for ACRI
## The 1st step of the derfinder pipeline is run as a snakemake pipeline. 

nohup snakemake --cores 8 --set-threads align=16 sort=16 index=16 &

### Also please keep config and snakefile files in the wdir. There should be a fastq folder with all fastq.gz files to be processed. IDs.txt should be a list of file ID prefixes, this file and libtype arguments are passed to config.yaml An environment running snakemake should be used to run. 
### COMBOseq fastqs are second-strand sequenced. 
### Illumina is FS sequences, therefore any smRNA pushed through exceRpt should be sense for comboseq and antisens for illumina
### This script uses a custom masked genome where sno, pi, miR and tRNA are considered contigs to avoid multimapping reads, other RNA types are mapped to genome and annotated with GENCODE.


## Followed by a helper script to merge counts into a RData object. 

Rscript Step2_derfinder_script_v3.R <inDir>

## Afterwards, the following command can be run to summarize results. 

R -e "rmarkdown::render('Step3_DE_edgeR.Rmd')" --args derfinder_out.Rda Design_template_1.txt 2 Responder Progression

#The arguments are 1 - the RData object from the helper script, 2 - The phenotype data table, 3 - The low-count filter (in CPM), 4 - the reference factor level, 5 - the treatment factor level.
#Using the counts from the excerpt_out folder, the supplied design matrix, the filter (CPM) and the reference and treatment levels, a HTML markdown is created.  
#Please refer to user eallain for any questions / details
