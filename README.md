> viGEN is a bioinformatics pipeline for the exploration of viral RNA in human NGS data. 
> In this tutorial, we provide and end to end workflow on how this pipeline can be used on an example data file.

## 1. Installation
- install RSEM 
```
git clone https://github.com/bli25broad/RSEM_tutorial.git
```
- install software (bowtie2 & rsem)
```
cd RSEM_tutorial/software
unzip bowtie2-2.2.6-source.zip
cd bowtie2-2.2.6
make -j 8
cd ..
tar -xzf RSEM-1.2.25.tar.gz
cd RSEM-1.2.25
make -j 8
make ebseq
cd ..
cd ..
#make -j 8 :8 threads. You can adjust this number according to your needs. 
#'make ebseq' compiles EBSeq and is only required if you want to use EBSeq for differential expression analysis.

```
#### Test RSEM successful installation.you can also add it to your $PATH.
```
$ ls
RSEM_tutorial  viGEN
$ ./RSEM_tutorial/software/RSEM-1.2.25/rsem-calculate-expression --version
Current version: RSEM v1.2.25
```   
  
    
      
## 2. Alignment to human reference （cost times）
preparation:   
1) fastq file(paired end):take SRR1946637 as an example in here.sorted in **input** directory . 
2) reference: download Rsem_ref_hg19.zip（https://drive.google.com/drive/folders/0B3-883ME4sP3dDF3ZllrN1JSWmM?usp=sharing） through google drive，unzip it,should be in the same directory with **input** directory .   
FOR more detail please see RSEM tutorial（https://github.com/bli25broad/RSEM_tutorial） in case you want to use your own ref file.

3) **bowtie**: there is an bowtie installing package followed but you can also use your own bowtie,type "whereis bowtie" and get it path,change it in command,then it also works.

```
$ whereis bowtie
bowtie: /path/to/bin/bowtie 
$ mkdir output
$ ./RSEM_tutorial/software/RSEM-1.2.25/rsem-calculate-expression --paired-end -p 20 --output-genome-bam --keep-intermediate files --bowtie-path /path/to/bin --append-names --estimate-rspd --time input/SRR1946637_1.fastq input/SRR1946637_2.fastq Rsem_ref_hg19/Rsem_ref_hg19 output/SRR1946637  
  
/path/to/bin/bowtie -q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 1000 -p 20 -a -m 200 -S Rsem_ref_hg19/Rsem_ref_hg19 -1 input/SRR1946637_1.fastq -2 input/SRR1946637_2.fastq | samtools view -S -b -o output/SRR1946637.temp/SRR1946637.bam -
Warning: Exhausted best-first chunk memory for read SRR1946637.1078913 HWI-ST1106:205:C0MLWACXX:4:1104:15798:26529 length=100/1 (patid 1078912); skipping read  
……  
  
```

option：  
  -p ：number of thread
  RSEM can used for paired end and single fastq,as well as SAM/BAM files, for more detail please check tutourial:http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html 
> SYNOPSIS  
     $rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name  
     $rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name  
     $rsem-calculate-expression [options] --sam/--bam [--paired-end] input reference_name sample_name  
     
Output files:   
In addition to the aligned BAM file (genome level and transcriptome level), this will generate the unaligned (unmapped) fastq files named SRR1946637_un_1.fq and SRR1946637_un_2.fq. They consist of the reads that did not align to the human reference.


## 3.Create viral reference
files download in https://drive.google.com/drive/folders/0B3-883ME4sP3Wm1FVjdVcEpfek0 from google drive.
- ***reference file***: writers obtained reference genomes of all known and sequenced human viruses obtained from NCBI (as of Sep 2015), and merged them into one file (referred to as the "viral reference file") in fasta file format. Merge all virus fasta file into one big fasta file called **viruses.fa**
- ***index file***: indexed the viral reference file, so that these files are ready for alignment tools Bowtie2 (folder name: virus.bowtie2.refB) or BWA (folder name: viral.bwa.ref)
> In case user is interested in creating a reference index the reference file on their own, this is the command to use: ```./bowtie2-build /Path/viruses.fa virus.bowtie2.refB```
- NCBI also allows to download information/annotation about these viruses from their web site. This information has been provided as Complete_Sequence_info.csv 

## 4.Align the unmapped fastq files to the viral reference  [bowtie2]

normal samtools process:  
```
```


## 5.Get genome level matrix file and find top viruses

