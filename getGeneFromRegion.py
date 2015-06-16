#!/usr/bin/env python

'''
Created on Aug 20, 2014
getGeneFromBlast
- The purpose of this script is to retrieve genes using a gff file 
using coordinates fromfrom a blast output file. 
- INPUT: 
    - blastoutput file
    - gff file
    - number of lines to stop
@author: mkatari
'''
#from Bio import SeqIO
from optparse import OptionParser

#import re
##################################################
# Do all the parsing of command line options here
##################################################
def parseArguments():
    parser= OptionParser()
    parser.add_option("-b", "--blastoutput", dest="blastoutput", default="test.region",
                      help="specify path to blast output file")
    parser.add_option("-g", "--gff", dest="gfffile", default="Mesculentav6.1.gene.gff3",
                      help="specify path to gff file")
    parser.add_option("-n", "--numhits", dest="numhits", default=1,
                      help="specify number of hits from the blast output file. If you want the best match then use 1")
    (options, args) = parser.parse_args()
    return options,parser

##################################################
# Do all the parsing of command line options here
##################################################
def readTabFile(infile, numhits):
    fh = open(infile, 'r')    
    
    lines = [i.strip().split('\t') for i in fh]    
    #print(lines)
    targets = {}

    for j in lines:
        if j[0].startswith('#'):
            pass
        else:

            if j[0] in targets.keys():
                if len(targets[j[0]]) == int(numhits):
                    continue
                else:
                    targets[j[0]].append([j[1], j[2], j[3]])

            else:
                targets[j[0]] = [[j[1], j[2], j[3]]]
            
    return targets

##################################################
# Do all the parsing of command line options here
##################################################
def getGenesFromGFF(infile):
    fh = open(infile, 'r')
    lines = [i.strip().split('\t') for i in fh]
    genes={}
    for j in lines:
        if j[0].startswith('#'):
            pass
        else:
            if j[2] != 'gene':
                pass
            else:
                if j[8] in genes.keys():
                    pass
                else:
                    genes[j[8]] = [j[0], j[3], j[4]]

    return genes

                

##################################################
##### MAIN #######################################
##################################################

if __name__ == '__main__':
    allarguments,parser = parseArguments()

    if allarguments.blastoutput is None and allarguments.gfffile is None and allarguments.numhits is None:
    parser.print_help()
    sys.exit(1)


    print "Loading Blastoutput files\n"
    targets = readTabFile(allarguments.blastoutput, allarguments.numhits)
    print "Done loading Blast\n"
    print "Loading GFF file\n";
    allgenes = getGenesFromGFF(allarguments.gfffile)
#    print(allgenes.keys())
#    print(targets)
    print "Done loading GFF\n";


    #foreach target
    for allk in targets.keys():
  #      print(allk)
        #each target has a list of lists. in case of top match, there should be only one list.
        for k in targets[allk]:
  #          print(k)
            targetchromosome = k[0].split("_")
            printtarget = (allk, k[0], k[1], k[2])
#            print "\t".join(printtarget)

            foundgene=0
            for g in allgenes.keys():
                printgene = (g.split(';')[1].split('=')[1],allgenes[g][0],allgenes[g][1],allgenes[g][2])
 #               print(printgene)
                if printtarget[1] == printgene[1]:
                    #print "chromosomes match"
                    if int(printtarget[2]) < int(printgene[3]):
                        #print "startcoor is less than geneendcoor"
                        if int(printtarget[3]) > int(printgene[2]):
                            #print "endcoor is greater than genestartcoor"
                            printtargetstring = [str(z) for z in printtarget]
                            printgenestring = [str(z) for z in printgene]
                            print "\t".join(printtargetstring + printgenestring)
                            foundgene = 1

            if foundgene == 0:
                nothingfound = ["nogenes", "0", "0", "0"]
                printtargetstring = [str(z) for z in printtarget]

                print "\t".join(printtargetstring + nothingfound)

    print "\n"
