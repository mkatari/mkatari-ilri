#!/usr/bin/env python

'''
Created on May 8, 2015
getFastaFromPosition
- The purpose of this script is to retrieve a fasta formatted file from a tab delimited file containing coordinate.
- Also this script will allow you to define how many basepairs surrounding the position should be retrieved.
- An alternate purpose is to define the start and stop of the sequence you want to retrieve.

- INPUT: 
    - fasta file with reference sequence (required)
    - position file with start and stop (if range) or just center
    - flanking sequence to print

- OUTPUT:
    - fasta file with sequence requested

@author: mkatari
'''
#from Bio import SeqIO

from optparse import OptionParser
import re

##################################################
# Do all the parsing of command line options here
##################################################

def parseArguments():
    parser = OptionParser()
    parser.add_option("-s", "--sequence", dest="seqfile",
                      help="specify path to fasta sequence file")
    parser.add_option("-p", "--positionfile", dest="position",
                      help="position file")
    parser.add_option("-f", "--flank", dest="flank",
                      help="amount of flanking sequence")
    (options, args) = parser.parse_args()
    return options

##################################################
# load fasta sequence
##################################################

def loadFastaSequence(seqfile):
    fh = open(seqfile, 'r')    
    seqs = {}
    sequence = []
    description = ""
    for i in fh.readlines():

        if i.startswith(">"):
            if len(sequence) != 0:
                totalsequence = "".join(sequence)
                seqs[description] = totalsequence

            sequence = []
            description = i.strip()
        else:
            sequence.append(i.strip())
    
    totalsequence = "".join(sequence)
    seqs[description] = totalsequence

    return seqs

##################################################
# Get fasta sequence
##################################################

def getSeqeunceFromPosition(allsequences, position, flank):
    fh = open(position, 'r')

    lines = [i.strip() for i in fh]
    #print(lines)

    for j in lines:
        eachrow = j.split(":")
        start = int(eachrow[1]) - flank
        end = int(eachrow[1]) + flank
        #print(eachrow, start, end)

        for k in allsequences.keys():
            tempk = "".join([">", k, " "])
            tempchr = "".join([">", eachrow[0], " "])
            matchfound = tempk.find(tempchr)
            #print(tempk, tempchr, matchfound)
            if tempk.find(tempchr) >= 0:
                newdisc = "".join([">", str(eachrow[0]), "_", str(start+1), "_", str(end)])
                print newdisc
                print allsequences[k][start:end]
                #print "\n"

        
##################################################
##### MAIN #######################################
##################################################

if __name__ == '__main__':
    allarguments = parseArguments()
    allsequences = loadFastaSequence(str(allarguments.seqfile))
    getSeqeunceFromPosition(allsequences, str(allarguments.position), int(allarguments.flank))
    #print(allsequences)

