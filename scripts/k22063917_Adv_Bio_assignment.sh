#!/bin/bash #

# Getting and installing Conda
# Download Miniconda3 latest version on Open Stack
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

#Then run this bash script to install miniconda3
bash Miniconda3-latest-Linux-x86_64.sh

# To reload our current shell session
source ~/.bashrc

#Then install the packages required

conda install fastqc
conda install trimmomatic
conda install bwa

conda install freebayes
conda install bedtools



conda install samtools
conda install picard

conda install vcflib

conda install snpEFF
conda install snpSift




#make the assignment directory just underneath the home
mkdir ~/k22063917_Adv_Bioinformatics_Assignment

#come back to home directory
cd

#change to the directory just made
cd k22063917_Adv_Bioinformatics_Assignment

#make new directories for the project
mkdir data meta results logs

#change directory to ~/k22063917_Adv_Bioinformatics_Assignment/data
cd data

#then make directories for untrimmed and trimmed fastq
mkdir untrimmed_fastq trimmed_fastq

#then change the directory to untrimmed_fastq
cd untrimmed_fastq

#then download the raw data
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

#decompress the binary raw fastsq files
bgzip -d NGS0001.R1.fastq.qz > NGS0001.R1.fastq
bgzip -d NGS0001.R2.fastq.qz > NGS0001.R2.fastq



#change directory to data (one above untrimmed_fastq)
cd ~/k22063917_Adv_Bioinformatics_Assignment/data

#then download the bed file
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

#perform quality assessment
cd ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq
fastqc -t 4 *.fastq

#make fastqc_untrimmed_reads directory then move the fastqc reports
mkdir ~/k22063917_Adv_Bioinformatics_Assignment/results/fastqc_untrimmed_reads
mv ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq/*fastqc* ~/k22063917_Adv_Bioinformatics_Assignment/results/fastqc_untrimmed_reads/

#perform trimmomatic
cd ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq



trimmomatic PE -threads 4 -phred33 ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq/NGS0001.R1.fastq ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq/NGS0001.R2.fastq -baseout ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10   TRAILING:25 MINLEN:50


#change directory to trimmed_fastq
cd ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/
#fastqc of paired trimmed fastq
fastqc -t 4 *P

# move fastqc reports of the trimmed reads to fastqc_trimmed_reads directory
cd ~/k22063917_Adv_Bioinformatics_Assignment/results/
mkdir fastqc_trimmed_reads
mv ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/*fastqc* ~/k22063917_Adv_Bioinformatics_Assignment/results/fastqc_trimmed_reads/


#now alignment

#download reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/k22063917_Adv_Bioinformatics_Assignment/data/

#make a new directory reference within data
# then move the reference file in therecd
cd ~/k22063917_Adv_Bioinformatics_Assignment/data/
mkdir reference
mv ~/k22063917_Adv_Bioinformatics_Assignment/data/hg19.fa.gz ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/


#index
bwa index ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa.gz

mkdir ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50  ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa.gz ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001.sam

# convert the sam file into bam format, sort it and generate an index using samtools

cd ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam
ls



#mark duplicates using picard tools

cd ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt


#using samtools, index the sorted & marked bam file, filter the reads according to the MAPQ quality(20) and filter on bitwise flag

cd ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data
samtools index NGS0001_sorted_marked.bam
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

#then index the filtered bam file again
samtools index NGS0001_sorted_filtered.bam

#then produce statistics flagstats, idxstats, depth of coverage
samtools flagstat NGS0001_sorted_filtered.bam
samtools idxstats NGS0001_sorted_filtered.bam
samtools stats NGS0001_sorted_filtered.bam

#get insert size
cd ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5
ls

#more statistics
cp ~/k22063917_Adv_Bioinformatics_Assignment/data/annotation.bed ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/annotation.bed
bedtools coverage -a ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001_sorted_filtered.bam -b ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/annotation.bed


#compresses some file to make space in the virtual machine
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq/NGS0001.R1.fastq
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/untrimmed_fastq/NGS0001.R2.fastq
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_1P
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_1U
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_2U
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/trimmed_fastq/NGS0001_trimmed_R_2P
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001.bam
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001.sam
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001_sorted.bam
bgzip ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001_sorted.bam.bai

# variant_Calling_filtering, using freebayes
#freebayes

zcat ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa.gz > ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa
samtools faidx ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa

freebayes --bam ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/k22063917_Adv_Bioinformatics_Assignment/data/reference/hg19.fa --vcf ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001.vcf

bgzip ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001.vcf
tabix -p vcf ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001.vcf.gz


#more statistics calculating depth of coverage
cp ~/k22063917_Adv_Bioinformatics_Assignment/data/annotation.bed ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/annotation.bed
bedtools coverage -a ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/NGS0001_sorted_filtered.bam -b ~/k22063917_Adv_Bioinformatics_Assignment/data/aligned_data/annotation.bed


#vcffilter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001.vcf.gz > ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered.vcfi


#filter the vcf file according to the annodation bed

bedtools intersect -header -wa -a ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered.vcf -b ~/k22063917_Adv_Bioinformatics_Assignment/data/annotation.bed > ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.vcf

bgzip ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.vcf

tabix -p vcf ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.vcf.gz


#unpack it and set annovar on Open Stack
cd ~
tar -zxvf annovar.latest.tar.gz

cd ~/annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp142 humandb/

cd ~/annovar
./convert2annovar.pl -format vcf4 ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.vcf.gz > ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.avinput

#generate csv output
cd ~/annovar
./table_annovar.pl ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.avinput humandb/ -buildver hg19 -out ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_annotated -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp142 -operation g,g,f,f,f,f -otherinfo -nastring . -csvout

#generate txt output
cd ~/annovar
./table_annovar.pl ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_hg19.avinput humandb/ -buildver hg19 -out ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_annotated -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro,avsnp142 -operation g,g,f,f,f,f -otherinfo -nastring .

#cut the txt output of annovar to leave columns that contain information desired, such as exonic or information from dbSNP
cut -f1,2,3,4,5,6,7,30 ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_annotated.hg19_multianno.txt > ~/k22063917_Adv_Bioinformatics_Assignment/results/NGS0001_filtered_annotated.hg19_multianno.txt.cut


#snpEff

cd ~/k22063917_Adv_Bioinformatics_Assignment/results
zcat NGS0001_filtered_hg19.vcf.gz > NGS0001_filtered_hg19.vcf

snpEff ann -download GRCh37.75 -canon -onlyProtein NGS0001_filtered_hg19.vcf >  NGS0001_filtered_hg19.ann.canon.vcf





