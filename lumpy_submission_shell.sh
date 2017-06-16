#!/bin/sh

# lumpy_variants.sh wrapper
# job array command, with taks 1-N
# When a single command in the array job is sent to a compute node,
# its task number is stored in the variable SGE_TASK_ID,
# so we can use the value of that variable to get the results we want:

# path to file with bam_id's
BAMs_to_PROCESS=/u/scratch2/c/colinpat/processed_BAMs_icnn.txt

# get the number of lines in txt file
number_bam=$(cat $BAMs_to_PROCESS | wc -l)

# set variables to pass to lumpy script
bam_directory=/u/nobackup/eeskin2/jhsul/bipolar/processed_BAMs_icnn
lumpy_folder=bipolar_lumpy
scratch=/u/scratch2/c/colinpat

# of jobs to process simultaneously 
JOBS=10


# submit jobs

qsub -t 1-$number_bam -tc $JOBS lumpy_variants.sh -f $BAMs_to_PROCESS -b $bam_directory -l $lumpy_folder -s $scratch
