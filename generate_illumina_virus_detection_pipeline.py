#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

# Example run:
# ./generate_illumina_virus_detection_pipeline.py -o /home/melody/Out -s /home/melody/ilm30.sa
# ./generate_illumina_virus_detection_pipeline.py -v -o /home/melody/testdir/vnnv_samples_test -s /home/melody/testdir/vnnv_test2.sa
# ./generate_illumina_virus_detection_pipeline.py -s /home/melody/testdir/test.sa

def main():
    user_directory = os.getcwd()
    home_directory = os.path.expanduser('~')

    # programs
    multiqc = home_directory + "/.local/bin/multiqc"
    fastqc  = "/usr/local/FastQC-0.11.9/fastqc"
    diamond = "/usr/local/diamond-2.0.11/diamond"
    seqtk   = "/usr/local/seqtk-1.3/seqtk"
    trinity = "/usr/local/trinityrnaseq-v2.13.1/Trinity"
    bwa     = "/usr/local/bwa-0.7.17/bwa"
    bowtie2 = "/usr/local/bowtie2-2.4.4/"
    epost   = "/usr/local/edirect-15.5/epost"
    efetch  = "/usr/local/edirect-15.5/efetch"
    blastn  = "/usr/local/ncbi-blast-2.12.0+/bin/blastn"
    generate_html = "/home/melody/testdir/generate_html.py"
    
    # databases
    uniprot90_virus_dmnd = "/usr/local/ref/uniprot/2021_02/virus/diamond/uniprot90_virus.dmnd"
    refseq_virus    = "/usr/local/ref/refseq/206/virus/blastdb/refseq.virus.fasta"

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_directory", help="Takes in path of directory that you want to create to store the output", type=str)
    parser.add_argument("-v", "--version", help="Tells you the version of script", action="store_true")
    parser.add_argument("-s", "--sample_file", help="Contains sample names and location of FASTQ files", type=str)
    args = parser.parse_args()

    if args.version:
        print("version 1.0")

    # Reading the file
    try:
        data_file = open(args.sample_file, "r")
        samples = data_file.readlines()
        data_file.close()
    except:
        print("You did not input a valid argument for -s, please see -h for more info.")
        return None

    if args.output_directory != None:
        if args.output_directory[0] != "/":
            args.output_directory = user_directory + "/" + args.output_directory
        print(f"Your directory is {args.output_directory}, your sample file is {args.sample_file}.")
    else:
        print(f"Default directory is QC_Out which has been created in your current directory, sample file is {args.sample_file}")
        args.output_directory = user_directory + "/QC_out"

    if args.version:
        print("version 1.0")

    # Create the directory
    if os.path.exists(args.output_directory):
        print(f"Warning: {args.output_directory} is an existing directory.")
    else:
        os.mkdir(args.output_directory)
        print(f"Directory {args.output_directory} created")

    # creates nested list for each sample
    sample_list = []
    for sample in samples:
        nested_sample_ls = []
        ls_samples = sample.split("\t")
        for sam in ls_samples:
            clean_sam = sam.rstrip("\n")
            nested_sample_ls.append(clean_sam)
        sample_list.append(nested_sample_ls)
      
    # creates a dictionary of samples
    d1 = {}
    for nested_sample_ls in sample_list:
        sample_name = nested_sample_ls[0]
        sample_R1 = nested_sample_ls[1]
        sample_R2 = nested_sample_ls[2]
        s = Sample(sample_name, sample_R1, sample_R2)
        d1[s.name] = s

    # change directory
    os.chdir(f"{args.output_directory}")
    cwd = os.getcwd()

    # creates directories for each sample
    for key in d1:
        sample_path = os.path.join(cwd, key)
        #if not os.path.exists(sample_path):
        Path(sample_path).mkdir(parents=True, exist_ok=True)
        print(f"Directory {key} created")

    # write makefile to user_directory
    os.chdir(f"{user_directory}")
    cwd = os.getcwd()
    print(f"Makefile will be created in {cwd}")

    WRITE_FLAG = "w"
    is_exist = os.path.isfile('./illumina_virus_detection_pipeline.mk')
    if (not is_exist):
        # Create a Makefile in the same location
        print("Created makefile")
        f = open(f"./illumina_virus_detection_pipeline.mk", WRITE_FLAG)
        f.close()
    else:
        print("Make file exists, overwritten")

    f = open(f"./illumina_virus_detection_pipeline.mk", WRITE_FLAG)

    f.write(".DELETE_ON_ERROR:\n\n")
    f.write("all:\t")

    for key in d1:
        f.write(f"{args.output_directory}/{key}/fastqc1.OK ")
        f.write(f"{args.output_directory}/{key}/diamond_fastq1.OK ")
        f.write(f"{args.output_directory}/{key}/subseq1.OK ")

        f.write(f"{args.output_directory}/{key}/fastqc2.OK ")
        f.write(f"{args.output_directory}/{key}/diamond_fastq2.OK ")
        f.write(f"{args.output_directory}/{key}/subseq2.OK ")

        f.write(f"{args.output_directory}/{key}/trinity_input.OK ")
        f.write(f"{args.output_directory}/{key}/trinity.OK ")
        f.write(f"{args.output_directory}/{key}/bwa.OK ")
        f.write(f"{args.output_directory}/{key}/refseqblast.OK ")
        f.write(f"{args.output_directory}/{key}/bowtie2.OK ")

        f.write(f"{args.output_directory}/{key}/html.OK ")

    f.write(f"{args.output_directory}/multiqc_output.OK")
    f.write("\n\n")

    # generates diamond reference database for indexing
#  f.write(f"\t{diamond} makedb --in {uniprot_database} -d {args.output_directory}/nr > {args.output_directory}/nr.log 2> {args.output_directory}/nr.err\n")
#  f.write(f"\ttouch {args.output_directory}/generate_db.OK\n\n")

    for key in d1:
        sample_object = d1[key]
        sampleR1_name = (((sample_object.fastq1.split("/"))[-1]).split("."))[0]
        sampleR2_name = (((sample_object.fastq2.split("/"))[-1]).split("."))[0]

        # run fastqc on fastq1
        f.write(f"{args.output_directory}/{key}/fastqc1.OK:\n")
        f.write(f"\t{fastqc} -o {args.output_directory}/{key} --extract {sample_object.fastq1} > {args.output_directory}/{key}/fastqc1.log 2> {args.output_directory}/{key}/fastqc1.err\n")
        f.write(f"\tmv {args.output_directory}/{key}/{sampleR1_name}_fastqc.html {args.output_directory}/{key}/fastqc1.html\n")
        f.write(f"\ttouch {args.output_directory}/{key}/fastqc1.OK\n\n")

        # run diamond on fastq1
        f.write(f"{args.output_directory}/{key}/diamond_fastq1.OK: \n")
        f.write(f"\tzcat {sample_object.fastq1} | seqtk seq -A > {args.output_directory}/{key}/input{key}1.fasta\n")
        f.write(f"\tcat {args.output_directory}/{key}/input{key}1.fasta | grep " +
                '">"' + f" | wc -l > {args.output_directory}/{key}/count1.txt\n")
        f.write(f"\t{diamond} blastx -d {uniprot90_virus_dmnd} -q {args.output_directory}/{key}/input{key}1.fasta -o {args.output_directory}/{key}/stitle{key}1.m8 -f 6 qseqid stitle > {args.output_directory}/{key}/diamond1.log 2> {args.output_directory}/{key}/diamond1.err\n")
        f.write(f"\tcat {args.output_directory}/{key}/stitle{key}1.m8 | cut -f1 | uniq -c | wc -l > {args.output_directory}/{key}/length_{key}1_diamond.txt\n")
        f.write(f"\ttouch {args.output_directory}/{key}/diamond_fastq1.OK\n\n")

        # subseq fastq1
        f.write(f"{args.output_directory}/{key}/subseq1.OK: {args.output_directory}/{key}/diamond_fastq1.OK\n")
        f.write(f"\t{seqtk} subseq {sample_object.fastq1} {args.output_directory}/{key}/stitle{key}1.m8 > {args.output_directory}/{key}/{key}_1.fq\n")
        f.write(f"\ttouch {args.output_directory}/{key}/subseq1.OK\n\n")

        # run fastqc on fastq2
        f.write(f"{args.output_directory}/{key}/fastqc2.OK:\n")
        f.write(f"\t{fastqc} -o {args.output_directory}/{key} --extract {sample_object.fastq2} > {args.output_directory}/{key}/fastqc2.log 2> {args.output_directory}/{key}/fastqc2.err\n")
        f.write(f"\tmv {args.output_directory}/{key}/{sampleR2_name}_fastqc.html {args.output_directory}/{key}/fastqc2.html\n")
        f.write(f"\ttouch {args.output_directory}/{key}/fastqc2.OK\n\n")

        # run diamond on fastq2
        f.write(f"{args.output_directory}/{key}/diamond_fastq2.OK: \n")
        f.write(f"\tzcat {sample_object.fastq2} | seqtk seq -A > {args.output_directory}/{key}/input{key}2.fasta\n")
        f.write(f"\tcat {args.output_directory}/{key}/input{key}2.fasta | grep " +
                '">"' + f" | wc -l > {args.output_directory}/{key}/count2.txt\n")
        f.write(f"\t{diamond} blastx -d {uniprot90_virus_dmnd} -q {args.output_directory}/{key}/input{key}2.fasta -o {args.output_directory}/{key}/stitle{key}2.m8 -f 6 qseqid stitle > {args.output_directory}/{key}/diamond2.log 2> {args.output_directory}/{key}/diamond2.err\n")
        f.write(f"\tcat {args.output_directory}/{key}/stitle{key}2.m8 | cut -f1 | uniq -c | wc -l > {args.output_directory}/{key}/length_{key}2_diamond.txt\n")
        f.write(f"\ttouch {args.output_directory}/{key}/diamond_fastq2.OK\n\n")

        # subseq fastq2
        f.write(f"{args.output_directory}/{key}/subseq2.OK: {args.output_directory}/{key}/diamond_fastq2.OK\n")
        f.write(f"\t{seqtk} subseq {sample_object.fastq2} {args.output_directory}/{key}/stitle{key}2.m8 > {args.output_directory}/{key}/{key}_2.fq\n")
        f.write(f"\ttouch {args.output_directory}/{key}/subseq2.OK\n\n")

        # generate trinity input files, subseq fastq1 and fastq2 again, using hits from both sides
        f.write(f"{args.output_directory}/{key}/trinity_input.OK: {args.output_directory}/{key}/subseq1.OK {args.output_directory}/{key}/subseq2.OK\n")
        f.write(f"\tcat {args.output_directory}/{key}/{key}_1.fq | grep " + '"@"' +
                f" | sed 's/ /\t/' | cut -f1 > {args.output_directory}/{key}/1_fq_accession.txt\n")
        f.write(f"\tcat {args.output_directory}/{key}/{key}_2.fq | grep " + '"@"' +
                f" | sed 's/ /\t/' | cut -f1 > {args.output_directory}/{key}/2_fq_accession.txt\n")
        f.write(f"\tcat {args.output_directory}/{key}/1_fq_accession.txt {args.output_directory}/{key}/2_fq_accession.txt > {args.output_directory}/{key}/combined_fq.txt\n")
        f.write(f"\tcat {args.output_directory}/{key}/combined_fq.txt | sort | uniq | sed -e 's/^@//' > {args.output_directory}/{key}/c.txt\n")
        f.write(f"\t{seqtk} subseq {sample_object.fastq1} {args.output_directory}/{key}/c.txt > {args.output_directory}/{key}/trinityR1_input.fq\n")
        f.write(f"\t{seqtk} subseq {sample_object.fastq2} {args.output_directory}/{key}/c.txt > {args.output_directory}/{key}/trinityR2_input.fq\n")
        f.write(f"\ttouch {args.output_directory}/{key}/trinity_input.OK\n\n")

        # trinity
        f.write(f"{args.output_directory}/{key}/trinity.OK: {args.output_directory}/{key}/trinity_input.OK\n")
        f.write(f"\t{trinity} --seqType fq --left {args.output_directory}/{key}/trinityR1_input.fq --right {args.output_directory}/{key}/trinityR2_input.fq --CPU 2 --max_memory 40G --output {args.output_directory}/{key}/trinity_assembly > {args.output_directory}/{key}/trinity_assembly.txt.log 2> {args.output_directory}/{key}/trinity_assembly.txt.err\n")
        f.write(f"\ttouch {args.output_directory}/{key}/trinity.OK\n\n")

        # mapping using bwa
        f.write(f"{args.output_directory}/{key}/bwa.OK: {args.output_directory}/{key}/trinity.OK\n")
        f.write(f"\t{bwa} index -a bwtsw {args.output_directory}/{key}/trinity_assembly.Trinity.fasta > {args.output_directory}/{key}/bwa_index.err 2> {args.output_directory}/{key}/bwa_index.log\n")
        f.write(f"\t{bwa} mem -t 2 -M {args.output_directory}/{key}/trinity_assembly.Trinity.fasta {args.output_directory}/{key}/trinityR1_input.fq {args.output_directory}/{key}/trinityR2_input.fq -o {args.output_directory}/{key}/out.sam> {args.output_directory}/{key}/bwa_mem.err 2> {args.output_directory}/{key}/bwa_mem.log\n")
        f.write(f"\tsamtools view -S -b {args.output_directory}/{key}/out.sam > {args.output_directory}/{key}/out.bam\n")
        f.write(f"\tsamtools sort {args.output_directory}/{key}/out.bam > {args.output_directory}/{key}/sorted.out.bam\n")
        f.write(f"\tsamtools coverage {args.output_directory}/{key}/sorted.out.bam > {args.output_directory}/{key}/coverage.txt\n")
        f.write(f"\ttouch {args.output_directory}/{key}/bwa.OK\n\n")

        # blastn against refseq
        # {args.output_directory}/{key}/trinity_assembly
        f.write(f"{args.output_directory}/{key}/refseqblast.OK: {args.output_directory}/{key}/bwa.OK\n")
        f.write(f"\t{blastn} -db {refseq_virus} -query {args.output_directory}/{key}/trinity_assembly.Trinity.fasta  -outfmt " +
                '"6 qacc stitle sacc pident"' + f" -out {args.output_directory}/{key}/{key}_blast.psl -max_target_seqs 1 > {args.output_directory}/{key}/blast_results{key}.txt.log 2> {args.output_directory}/{key}/blast_results{key}.txt.err\n")
        f.write(f"\tcut -f2 {args.output_directory}/{key}/{key}_blast.psl | sort | uniq -c > {args.output_directory}/{key}/blast_{key}.txt\n")
        f.write(f"\ttouch {args.output_directory}/{key}/refseqblast.OK\n\n")

        # map reads to hits using bowtie2
        f.write(f"{args.output_directory}/{key}/bowtie2.OK: {args.output_directory}/{key}/refseqblast.OK\n")
        f.write(f"\tcat {args.output_directory}/{key}/{key}_blast.psl | cut -f3 | sort | uniq -c > {args.output_directory}/{key}/bowtie_accession.txt\n")
        f.write(f"\t{epost} -db nucleotide -input {args.output_directory}/{key}/bowtie_accession.txt | {efetch} -format fasta > {args.output_directory}/{key}/references.fna\n")
        f.write(f"\t{bowtie2}bowtie2-build {args.output_directory}/{key}/references.fna {args.output_directory}/{key}/bowtie_database > {args.output_directory}/{key}/bowtie2_index.txt.log 2> {args.output_directory}/{key}/bowtie2_index.txt.err\n")
        f.write(f"\t{seqtk} seq -F '#' {args.output_directory}/{key}/trinity_assembly.Trinity.fasta > {args.output_directory}/{key}/trinity_assembly.Trinity.fastq\n")
        f.write(f"\t{bowtie2}bowtie2 -x  {args.output_directory}/{key}/bowtie_database -U {sample_object.fastq1},{sample_object.fastq2} -S {args.output_directory}/{key}/bowtie2_R1.sam > {args.output_directory}/{key}/bowtie2.txt.log 2> {args.output_directory}/{key}/bowtie2.txt.err\n")
        f.write(f"\ttouch {args.output_directory}/{key}/bowtie2.OK\n\n")

        # generate html and js
        # ./generate_html.py -o /home/melody/testdir/QC_out/6_BTV-SISPA/html -v -q /home/melody/testdir/QC_out/6_BTV-SISPA/
        f.write(f"{args.output_directory}/{key}/html.OK: {args.output_directory}/{key}/refseqblast.OK\n")
        f.write(f"\t{generate_html} -o {args.output_directory}/{key}/html -v -q {args.output_directory}/{key}/ > {args.output_directory}/{key}/html.txt.log 2> {args.output_directory}/{key}/html.txt.err\n")
        f.write(f"\ttouch {args.output_directory}/{key}/html.OK\n\n")

    # multiqc
    f.write(f"{args.output_directory}/multiqc_output.OK:")
    for key in d1:
        f.write(f" {args.output_directory}/{key}/diamond_fastq1.OK {args.output_directory}/{key}/diamond_fastq2.OK")
    f.write("\n")
    f.write(f"\t{multiqc} {args.output_directory} -o {args.output_directory}/multiqc_output/\n")
    f.write(f"\ttouch {args.output_directory}/multiqc_output.OK\n\n")

    # clean
    f.write("clean:\n")
    f.write(f"\trm -f {args.output_directory}/nr* ")
    for key in d1:
        f.write(f"{args.output_directory}/{key}/*.log {args.output_directory}/{key}/*.err {args.output_directory}/{key}/input{key}1.fasta {args.output_directory}/{key}/input{key}2.fasta {args.output_directory}/{key}/1_fq_accession.txt {args.output_directory}/{key}/2_fq_accession.txt {args.output_directory}/{key}/combined_fq.txt {args.output_directory}/{key}/trinityR1_input.fq {args.output_directory}/{key}/trinityR2_input.fq {args.output_directory}/{key}/bowtie_accession.txt {args.output_directory}/{key}/references.fna {args.output_directory}/{key}/bowtie_database* {args.output_directory}/{key}/trinity_assembly.Trinity.fastq ")
    f.close()

# creates sample class
class Sample(object):
    def __init__(self, name, fastq1, fastq2):
        self.name = name
        self.fastq1 = fastq1
        self.fastq2 = fastq2

main()