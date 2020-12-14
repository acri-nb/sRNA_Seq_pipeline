# sRNA_Seq_pipeline


The pipeline work in three steps:

Step1 : Fastq processing, with adapter trimming plus mapping

Step2: Derfinder pipeline to get the count matrix and annotation of chromosome locations

Step3: Differential expression steps.



Both Step2 and Step3 happens in R, however the step2 is important to be run in a server, since it needs processing capabilities.
