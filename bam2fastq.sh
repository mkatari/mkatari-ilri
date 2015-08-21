#!/bin/env bash
#SBATCH -p batch
#SBATCH -J sam2fastq
#SBATCH -n 1
#SBATCH -w taurus

mkdir /var/scratch/$USER

TMPDIR=`mktemp -d --tmpdir=/var/scratch/$USER/`
PWD=`pwd`

INPUT=$1

#Name is assumed to be the first word before any period
NAME=`ls $INPUT | cut -f 1 -d '.'`
LEFT=${NAME}_1.fastq
RIGHT=${NAME}_2.fastq

cp $INPUT ${TMPDIR}/

echo "Ready to sam2fastq"

cd $TMPDIR

java -Djava.io.tmpdir=/var/scratch/mkatari/tmp -jar /export/apps/picard-tools/1.112/SamToFastq.jar \
    INPUT=${TMPDIR}/${INPUT} \
    FASTQ=${TMPDIR}/${LEFT} \
    SECOND_END_FASTQ=${TMPDIR}/${RIGHT} \

mv ${TMPDIR}/${LEFT} ${PWD}/${LEFT}

mv ${TMPDIR}/${RIGHT} ${PWD}/${RIGHT}


cd $PWD

echo "Cleaning up"
rm -rf ${TMPDIR}


echo "Done sam2fastq"


