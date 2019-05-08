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
  
    
      
## 2-1. Alignment to human reference （cost times）
preparation:   
1) fastq file(paired end):take SRR1946637 as an example in here.sorted in **input** directory . 
2) reference: download Rsem_ref_hg19.zip（https://drive.google.com/drive/folders/0B3-883ME4sP3dDF3ZllrN1JSWmM?usp=sharing） through google drive，unzip it,should be in the same directory with **input** directory .   
FOR more detail please see RSEM tutorial（https://github.com/bli25broad/RSEM_tutorial） in case you want to use your own ref file.(rsem-prepare-reference + ref file + gtt file)

3) **bowtie**: there is an bowtie installing package followed but you can also use your own bowtie,type "whereis bowtie" and get it path,change it in command,then it also works.

```
$ whereis bowtie
bowtie: /path/to/bin/bowtie 
$ mkdir output
$ ./RSEM_tutorial/software/RSEM-1.2.25/rsem-calculate-expression --paired-end -p 20 --output-genome-bam --keep-intermediate-files --bowtie-path /path/to/bin --append-names --estimate-rspd --time input/SRR1946637_1.fastq input/SRR1946637_2.fastq Rsem_ref_hg19/Rsem_ref_hg19 output/SRR1946637  
  
/path/to/bin/bowtie -q --phred33-quals -n 2 -e 99999999 -l 25 -I 1 -X 1000 -p 20 -a -m 200 -S Rsem_ref_hg19/Rsem_ref_hg19 -1 input/SRR1946637_1.fastq -2 input/SRR1946637_2.fastq | samtools view -S -b -o output/SRR1946637.temp/SRR1946637.bam -
Warning: Exhausted best-first chunk memory for read SRR1946637.1078913 HWI-ST1106:205:C0MLWACXX:4:1104:15798:26529 length=100/1 (patid 1078912); skipping read  
……  

#single end
#./RSEM_tutorial/software/RSEM-1.2.25/rsem-calculate-expression -p 20 --output-genome-bam --keep-intermediate-files --bowtie-path /path/to/bin --append-names --estimate-rspd --time input/SRR5275286.fastq Rsem_ref_hg19/Rsem_ref_hg19 output/SRR5275286
  
```

option：  
  -p ：number of thread  
  --keep-intermediate-files : continue after unexcepted stop  
  RSEM can used for paired end and single fastq,as well as SAM/BAM files, for more detail please check tutourial:http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html 

     
Output files:   
In addition to the aligned BAM file (genome level and transcriptome level), this will generate the unaligned (unmapped) fastq files named SRR1946637_un_1.fq and SRR1946637_un_2.fq. They consist of the reads that did not align to the human reference.

  
if error:   
> Invalid number of arguments!  
NAME  
    rsem-calculate-expression   
    
It reminds you that your something wrong with your arguments or option,please ensure you have set right argument and carefully check all -option spelling  

> SYNOPSIS  
     $rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name  
     $rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name  
     $rsem-calculate-expression [options] --sam/--bam [--paired-end] input reference_name sample_name  


## 2-2-1.Create viral reference
files download in https://drive.google.com/drive/folders/0B3-883ME4sP3Wm1FVjdVcEpfek0 from google drive.
- ***reference file***: writers obtained reference genomes of all known and sequenced human viruses obtained from NCBI (as of Sep 2015), and merged them into one file (referred to as the "viral reference file") in fasta file format. Merge all virus fasta file into one big fasta file called **viruses.fa**
- ***index file***: indexed the viral reference file, so that these files are ready for alignment tools Bowtie2 (folder name: virus.bowtie2.refB) or BWA (folder name: viral.bwa.ref)
> In case user is interested in creating a reference index the reference file on their own, this is the command to use: ```./bowtie2-build /Path/viruses.fa virus.bowtie2.refB```
- NCBI also allows to download information/annotation about these viruses from their web site. This information has been provided as Complete_Sequence_info.csv 

## 2-2-2.Align the unmapped fastq files to the viral reference 
- ***software***: [bowtie2] ,[Samtools]
- ***process***: normal samtools process, you can use the shell script provided ```viral.pipeline_public_final.sh``` that encompasses all of these steps.

```
$ ls
3_Get.viral.reference  bowtie.sh    count.sh   nohup.out  output RSEM_tutorial file.list  normal.txt  Rsem_ref_hg19  viGEN  viral.pipeline_public_final.sh
$ ls Rsem_ref_hg19/
Rsem_ref_hg19.1.ebwt  Rsem_ref_hg19.4.ebwt   Rsem_ref_hg19.idx.fa      Rsem_ref_hg19.seq             RSEMRef.stderr.txt.docx
Rsem_ref_hg19.2.ebwt  Rsem_ref_hg19.chrlist  Rsem_ref_hg19.rev.1.ebwt  Rsem_ref_hg19.ti              RSEMRef.stdout.txt.docx
Rsem_ref_hg19.3.ebwt  Rsem_ref_hg19.grp      Rsem_ref_hg19.rev.2.ebwt  Rsem_ref_hg19.transcripts.fa
$ ls 3_Get.viral.reference/
viGen_ref.txt  viral.bowtie2.ref  viral.bwa.ref  viruses.dict  viruses.fa  viruses.fa.fai
$head file.list
SRR1946637  
SRR1946638  
SRR1946639  
```

error：Segmentation fault (core dumped) (ERR): bowtie2-align exited with value 139  
it happens when you deal with a group project,did not effect sigle task,you need to update your samtools version to fix this problem.  
error：Segmentation fault (core dumped) (ERR): bowtie2-align exited with value 1
check your option and parameter please

***you can choose method 2-1(align to humman ref first and get unmap.fastq, then re-align to viral reference to get an bam file )  or method 2-2 (align to vrial refernece directory to get the bam file),or both***

## change an reference file

```
#build index
./bowtie2-build /Path/to/new_viruses.fa new_virus.bowtie2.ref
Wrote 12362996 bytes to primary EBWT file: hervquant_virus.bowtie2.ref.rev.1.bt2
Wrote 5969992 bytes to secondary EBWT file: hervquant_virus.bowtie2.ref.rev.2.bt2
Re-opening _in1 and _in2 as input streams
Returning from Ebwt constructor
Headers:
    len: 23879943
    bwtLen: 23879944
    sz: 5969986
    bwtSz: 5969986
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 20
    eftabSz: 80
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 1492497
    offsSz: 5969988
    lineSz: 64
    sideSz: 64
    sideBwtSz: 48
    sideBwtLen: 192
    numSides: 124375
    numLines: 124375
    ebwtTotLen: 7960000
    ebwtTotSz: 7960000
    color: 0
    reverse: 1
Total time for backward call to driver() for mirror index: 00:00:12
```

## output and statistic analysis

output file
```
ll output/SRR1946637/
total 44
log.bowtie2.txt
SRR1946637.bowtie2.idxstats.txt
system.err.bowtie2.txt
$ head output/SRR1946637/*.id*

##header:a virus sequence name, sequence length, # mapped reads and # unmapped reads
gi|10313991|ref|NC_002549.1|	18959	0	0
gi|106060735|ref|NC_001959.2|	7654	32 0
gi|10937870|ref|NC_001796.2|	15462	0	0
gi|109390382|ref|NC_008188.1|	7263	0	0
gi|109390389|ref|NC_008189.1|	7259	450	0
gi|110645916|ref|NC_001401.2|	4679	2	0
gi|11528013|ref|NC_001563.2|	10962	0	0
gi|11545722|ref|NC_002617.1|	15186	0	
gi|119952252|ref|NC_008719.1|	10793	3	0
gi|119952254|ref|NC_008718.1|	10510	0	0
```
use ```$bash count.sh``` can get the sum of column 3(#mapped read) of *.bowtie2.idxstats.txt* file

