> viGEN is a bioinformatics pipeline for the exploration of viral RNA in human NGS data. 
> In this tutorial, we provide and end to end workflow on how this pipeline can be used on an example data file.

## 1. Installation
- install RSEM 
```
git clone https://github.com/bli25broad/RSEM_tutorial.git
```
- install software (bowtie2 & rsem
```
cd software
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
*** Test RSEM successful installation.you can also add it to your $PATH.
```
$ ls
RSEM_tutorial  viGEN
$ ./RSEM_tutorial/software/RSEM-1.2.25/rsem-calculate-expression --version
Current version: RSEM v1.2.25
``` 
  
