#!/bin/bash
#$ -pe shared 6
#$ -l h_rt=20:00:00
#$ -l h_data=2G
#$ -o /u/home/c/colinpat/batch_scripts/lumpy_batch/lumpy_output/
#$ -e /u/home/c/colinpat/batch_scripts/lumpy_batch/lumpy_output/

while getopts f:b:l:s:v: option
do
 case "${option}"
 in
 f) BAMs_to_PROCESS=${OPTARG};;
 b) bam_directory=${OPTARG};;
 l) lumpy_folder=${OPTARG};;
 s) scratch=${OPTARG};;
 v) lumpy_vcf=$OPTARG;;
 esac
done

mkdir -p ${scratch}/${lumpy_folder}
cd ${scratch}/${lumpy_folder}
. /u/local/Modules/default/init/modules.sh
module load python/2.7.3
module load samtools/1.3.1
module load bedtools
module load R

sample_name=$(cat $BAMs_to_PROCESS | head -${SGE_TASK_ID} | tail -1 )

echo $SGE_TASK_ID 
echo $BAMs_to_PROCESS
echo $sample_name
echo $bam_directory
echo $lumpy_folder
echo $scratch

/u/home/c/colinpat/Lumpy/svtyper-0.1.4/svtyper -B ${bam_directory}/${sample_name}/${sample_name}.dedup.realigned.recal.bam -i ${scratch}/${lumpy_folder}/${lumpy_vcf} --max_reads 1000 > ${scratch}/${lumpy_folder}/${sample_name}.svtyper.vcf
