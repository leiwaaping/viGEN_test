#!/bin/bash

# Align to viral reference using Bowtie2
functionBowtie2() {
    echo "Start - Align to viral reference using Bowtie2 local algin"
    bowtie2 -x ${BOWTIE2_REF} --local -1 ${I_FOLDER}${ID}_1.fastq -2 ${I_FOLDER}${ID}_2.fastq --threads ${N_THREAD} --al ${O_FOLDER}${ID}/${ID}.bowtie2.align --un ${O_FOLDER}${ID}/${ID}.Bowtie2.unaligned -S ${O_FOLDER}${ID}/${ID}.bowtie2.sam --sensitive
    echo "End "
}

#Convert SAM file to BAM file
functionConvertSamtoBam() {
    echo "Start - convert sam to bam file using Samtools" 
#    samtools view -bS -@ $N_THREAD ${O_FOLDER}${ID}/${ID}.bowtie2.sam > ${O_FOLDER}${ID}/${ID}.bowtie2.bam &&
    samtools view -bS -@ $N_THREAD ~/../../data/HPli_data/${ID}.bowtie2.sam > ${O_FOLDER}${ID}/${ID}.bowtie2.bam &&
    echo "End " 
    #save space
    rm ${O_FOLDER}${ID}/${ID}.bowtie2.sam
}

#Sort BAM coordinate wise
functionSortBam() {
    echo "Start - Samtools sort bam file" 
    samtools sort -@ $N_THREAD ${O_FOLDER}${ID}/${ID}.bowtie2.bam -o ${O_FOLDER}${ID}/${ID}.bowtie2.sorted.bam &&
    echo "End " 
    rm ${O_FOLDER}${ID}/${ID}.bowtie2.bam 
}

#Index BAM file
functionIndexBam() {
    echo "Start - Samtools index bam file"
    samtools index -@ $N_THREAD ${O_FOLDER}${ID}/${ID}.bowtie2.sorted.bam &&
    echo "End "
}

# Samtools idx -  genome level counts
functionIdxStatsBam() {
    echo "Start - Samtools Idx"
    samtools idxstats -@ $N_THREAD ${O_FOLDER}${ID}/${ID}.bowtie2.sorted.bam > ${O_FOLDER}${ID}/${ID}.bowtie2.idxstats.txt &&
    echo "End "
 }

########################################
# Values to set by the user
N_THREAD=20
I_FOLDER="~/../../data/HPli_data/SRA/"
O_FOLDER="output/"

LOG=$O_FOLDER$ID/"pipeline_log.txt"
ERR2="system.err.bowtie2.txt"
ERR3="log.bowtie2.txt"
#SAMTOOLS="~/.conda/envs/bioinenv/bin/"
#BOWTIE2_SW="~/.conda/envs/bioinenv/bin/"
BOWTIE2_REF="./3_Get.viral.reference/viral.bowtie2.ref/virus.bowtie2.ref"
		
for x in $(cat samplenames.txt)
do
    echo $x
    ID=${x%*}

    echo "-----------------------------" 
    echo "precessing :$ID" 

    START=$(date +%s) 
    
    mkdir $O_FOLDER$ID

    functionBowtie2 &&

    functionConvertSamtoBam &&

    functionSortBam &&

    functionIndexBam &&

    functionIdxStatsBam &&

    #functionFlagstatBam &&
    
    END=$(date +%s)
    DIFF=$((($END-$START)/60))
    echo "It took $DIFF minutes"
done
