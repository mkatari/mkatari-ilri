#!/bin/env bash
#SBATCH -p batch
#SBATCH -J bowtie2
#SBATCH -n 4
#SBATCH -w taurus

#Creating a temp directory.
#The script should be executed from the folder that contains files.

mkdir /var/scratch/$USER

TMPDIR=`mktemp -d --tmpdir=/var/scratch/$USER/`
PWD=`pwd`

#getting and setting variables
DATABASE=$1
LEFT=$2
RIGHT=$3
OUTPUT=$4
SAM=${OUTPUT}.sam
BAM=${OUTPUT}.bam
SORTED=${OUTPUT}_sorted.bam

#move to temp directory
cp ${LEFT} ${TMPDIR}/
cp ${RIGHT} ${TMPDIR}/

echo "Ready to align using bowtie2"

#run bowtie

cd $TMPDIR

module load bowtie2/2.2.3

bowtie2 --local -x ${DATABASE} \
-p 4 \
-1 ${TMPDIR}/${LEFT} \
-2 ${TMPDIR}/${RIGHT} \
-S ${TMPDIR}/${SAM} 

echo "Converting to bam"

#convert to bam
module load samtools/0.1.19

samtools view -bS ${TMPDIR}/$OUTPUT > ${TMPDIR}/${BAM}


echo "Sorting bam"
#sort using picard
java –Djava.io.tmpdir=$TMPDIR \
–jar /export/apps/picard-tools/1.112/SortSam.jar \
INPUT=${TMPDIR}/${BAM} \
OUTPUT=${TMPDIR}/${SORTED} \
SORT_ORDER=coordinate

echo "Moving results to pwd."
mv ${TMPDIR}/${SORTED} ${PWD}/${SORTED}

cd $PWD

echo "Cleaning up"
rm -rf ${TMPDIR}

echo "Done with bowtie2"


