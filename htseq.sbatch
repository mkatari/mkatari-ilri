#!/bin/env bash
#SBATCH -p batch
#SBATCH -J htseqcount
#SBATCH -n 1
#SBATCH -w taurus

#Creating a temp directory.
#The script should be executed from the folder that contains files.

mkdir /var/scratch/$USER
TMPDIR=`mktemp -d --tmpdir=/var/scratch/$USER/`
PWD=`pwd`

#getting and setting variables
BAMFILE=$1
GFF=$2
OUTPUT=$3

SORTED=${BAMFILE}.sorted
SORTEDBAM=${SORTED}.bam
#SAMFILE=${BAMFILE}.sam

#move to temp directory
cp ${BAMFILE} ${TMPDIR}/
cp ${GFF} ${TMPDIR}/

echo "Ready to count using htseqcount"

#run samtools

cd $TMPDIR

module load samtools/0.1.19
samtools sort -n ${TMPDIR}/${BAMFILE} ${TMPDIR}/${SORTED}

ls ${TMPDIR}

#samtools view -h ${TMPDIR}/${SORTED}.bam -o ${TMPDIR}/${SAMFILE}

#run htseq
#module load htseq/0.5.4p5

module load htseq/0.6.1

htseq-count -f bam -t gene -i Name -m union -r name -s no ${TMPDIR}/${SORTEDBAM} $GFF > ${TMPDIR}/${OUTPUT}

echo "Moving results to pwd."
mv ${TMPDIR}/${OUTPUT} ${PWD}/${OUTPUT}

cd $PWD

echo "Cleaning up"
rm -rf ${TMPDIR}

echo "Done with htseq count"

