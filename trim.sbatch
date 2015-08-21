#!/bin/env bash
#SBATCH -p batch
#SBATCH -J trimmomatic
#SBATCH -n 1
#SBATCH -w taurus

LEFT=$1
RIGHT=$2

TMPDIR=`mktemp -d --tmpdir=/var/scratch/mkatari/`
PWD=`pwd`

cp $LEFT ${TMPDIR}/
cp $RIGHT ${TMPDIR}/

OUTPUT1PAIR=${LEFT}_pair.fastq
OUTPUT1UNPAIR=${LEFT}_unpair.fastq
OUTPUT2PAIR=${RIGHT}_pair.fastq
OUTPUT2UNPAIR=${RIGHT}_unpair.fastq

echo "Ready to trim"

java -jar /home/mkatari/apps/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 ${TMPDIR}/$LEFT \
	${TMPDIR}/$RIGHT ${TMPDIR}/$OUTPUT1PAIR ${TMPDIR}/$OUTPUT1UNPAIR ${TMPDIR}/$OUTPUT2PAIR ${TMPDIR}/$OUTPUT2UNPAIR \
	HEADCROP:15 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


mv ${TMPDIR}/${OUTPUT1PAIR} ${PWD}/${OUTPUT1PAIR}
mv ${TMPDIR}/${OUTPUT1UNPAIR} ${PWD}/${OUTPUT1UNPAIR}
mv ${TMPDIR}/${OUTPUT2PAIR} ${PWD}/${OUTPUT2PAIR}
mv ${TMPDIR}/${OUTPUT2UNPAIR} ${PWD}/${OUTPUT2UNPAIR}


echo "Done trimming"


