#!/bin/bash
#$ -pe shared 4
#$ -l h_rt=10:00:00
#$ -l h_data=2G
#$ -o /u/home/c/colinpat/batch_scripts/lumpy_batch/lumpy_output/
#$ -e /u/home/c/colinpat/batch_scripts/lumpy_batch/lumpy_output/

while getopts f:b:l:s: option
do
 case "${option}"
 in
 f) BAMs_to_PROCESS=${OPTARG};;
 b) bam_directory=${OPTARG};;
 l) lumpy_folder=${OPTARG};;
 s) scratch=$OPTARG;;
 esac
done

echo $SGE_TASK_ID 

sample_name=$(cat $BAMs_to_PROCESS | head -${SGE_TASK_ID} | tail -1 )

echo $BAMs_to_PROCESS
echo $sample_name
echo $bam_directory
echo $lumpy_folder
echo $scratch

mkdir -p ${scratch}/${lumpy_folder}
cd ${scratch}/${lumpy_folder}
. /u/local/Modules/default/init/modules.sh
module load python/2.7.13
module load samtools/1.3.1
module load bedtools
module load R


/u/home/c/colinpat/Lumpy/extract_sv_reads-1.1.2/build/bin/extract-sv-reads -e -i ${bam_directory}/${sample_name}/${sample_name}.dedup.realigned.recal.bam -s ${scratch}/${lumpy_folder}/${sample_name}.splitters.bam -d ${scratch}/${lumpy_folder}/${sample_name}.discord.bam
samtools index ${scratch}/${lumpy_folder}/${sample_name}.splitters.bam
samtools index ${scratch}/${lumpy_folder}/${sample_name}.discord.bam
stats=$(samtools view -r ${sample_name} ${bam_directory}/${sample_name}/${sample_name}.dedup.realigned.recal.bam | tail -n+100000 | /u/home/c/colinpat/Lumpy/lumpy-sv/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o ${sample_name}.histo | tail -n 1 | sed 's/[[:space:]]/:/' | awk '{split($0,stats,":"); print stats[2],stats[4]}')
MEAN=$(echo $stats | cut -f1 -d ' ')
STDEV=$(echo $stats | cut -f2 -d ' ')
/u/home/c/colinpat/Lumpy/lumpy-sv-v0.2.13/bin/lumpy -P -mw 4 -tt 0 -pe id:${sample_name},bam_file:${sample_name}.discord.bam,histo_file:${sample_name}.histo,mean:${MEAN},stdev:${STDEV},read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:${sample_name},bam_file:${sample_name}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > ${sample_name}.lumpy.vcf
