configfile: "config.yaml"
#include: "utils.py"
#Command line syntax:
#snakemake -s snakefile.smk --cores 8 --config prefix= lib=
class MyException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

with open(config['IDs'], 'r') as file:
    # Read all the lines of the file into a list
    lines = file.read().splitlines()
rule all: 
    input:
        expand("fastq/{prefix}.fastq.gz", prefix=lines),
        expand("trimfastq/{prefix}/{prefix}_trimmed.fq.gz", prefix = lines),
        expand("trimfastq/{prefix}", prefix = lines),
        expand("reportfalco/{prefix}", prefix=lines),
        expand("sam/{prefix}Aligned.out.bam", prefix = lines),
        expand("sam/{prefix}.s.bam", prefix=lines),
        expand("sam/{prefix}.s.bam.bai", prefix=lines) 

#Checking fastq quality
rule testquatility:
    input:
        fq = "fastq/{prefix}.fastq.gz"
    output:
        falrep = directory("reportfalco/{prefix}")
    shell:
        """
        /home/acri/tools/Falco/falco/bin/falco {input.fq} -o {output.falrep}
        """
#Trimming
rule trimming:
    input:
        fq_in = "fastq/{prefix}.fastq.gz"
    output:
        trimdir = directory("trimfastq/{prefix}"),
        fq_out = "trimfastq/{prefix}/{prefix}_trimmed.fq.gz"
    run:
        if config['Lib']=='nextflex':
            shell("/home/acri/tools/TrimGalore-0.6.6/trim_galore -a TGGAATTCTCGGGTGCCAAGG --clip_R1 4 --three_prime_clip_R1 4 --length 18 {input.fq_in} --output_dir {output.trimdir}")
        elif config['Lib']=='comboseq':
            shell("/home/acri/tools/TrimGalore-0.6.6/trim_galore -a AAAAAAAA --length 15 {input.fq_in} --output_dir {output.trimdir}")
        elif config['Lib']=='trilink':
            shell("/home/acri/tools/TrimGalore-0.6.6/trim_galore -a TCACCGACTGCCATAGAG --small_rna {input.fq_in} --output_dir {output.trimdir}")
        else:
            raise MyException("Library strategy is expected to be either: trilink, comboseq or nextflex; but is" + str(config['Lib']))
            with MyException as e:
                print("\033[91mError when we try to execute trimmimg :", e.message, "\033[0m")

#Alignment
rule align:
    input:
        ref="/references/DERFINDER_ref/STAR_INDEX",
        trim_in="trimfastq/{prefix}/{prefix}_trimmed.fq.gz"
    params:
        pr="sam/{prefix}"
    output:
        "sam/{prefix}Aligned.out.bam"
    run:
        if config['Lib']=='nextflex' or config['Lib']=='trilink':
            shell("/home/acri/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir {input.ref} --readFilesIn {input.trim_in} --readFilesCommand gunzip -c --outFileNamePrefix {params.pr} --outFilterScoreMinOverLread 0 --outFilterMatchNmin 16 --outFilterMatchNminOverLread 0 --outFilterMismatchNoverLmax 0.05 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMtype BAM Unsorted")
        elif config['Lib']=='comboseq':
            shell("/home/acri/tools/STAR-2.7.9a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir {input.ref} --readFilesIn {input.trim_in} --readFilesCommand gunzip -c --outFileNamePrefix {params.pr} --outFilterScoreMinOverLread 0 --outFilterMatchNmin 15 --outFilterMatchNminOverLread 0.9 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0.3 --alignIntronMax 1 --outSAMtype BAM Unsorted")
        else:
            raise MyException("Library strategy is expected to be either: trilink, comboseq or nextflex; but is" + str(config['Lib']))
            with MyException as e:
                print("\033[91mError when we try to execute trimmimg :", e.message, "\033[0m")


#Housekeeping
rule sort:
    input:
        al_in="sam/{prefix}Aligned.out.bam",
    output:
        sambam="sam/{prefix}.s.bam"
    shell:
        """
        samtools sort -@ 16 -o {output.sambam} {input.al_in}; samtools index {output.sambam}; rm {input.al_in}
        """

#index reads
rule index:
    input:
        al_in="sam/{prefix}.s.bam",
    output:
        sambam="sam/{prefix}.s.bam.bai"
    shell:
        """
        samtools index {input.al_in}
        """

