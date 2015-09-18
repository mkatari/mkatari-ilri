#!/usr/bin/env python

'''
Created on June 6, 2015
goHyperG.py
- The purpose of this script is to calculate GO enrichment using hypergeometric test

- INPUT:
    - genelist a file with a list of genes (required)
    - a goterm association file (it doesn't have to be go-terms
    - a file with description of the go-terms.

- OUTPUT:
    - a tab delimited file / table summarizing the terms and their p-values

@author: mkatari
'''
__author__ = 'manpreetkatari'

from optparse import OptionParser
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import sys

##################################################
# Do all the parsing of command line options here
##################################################

def parseArguments():
    parser = OptionParser()
    parser.add_option("-g", "--genelist", dest="genelistfile", default="cassava_genelist.txt",
                      help="specify path to file with list of genes")
    parser.add_option("-s", "--species", dest="species", default="cassavaV6",
                      help="specify the species to use")
    parser.add_option("-t", "--term", dest="term", default="ALL",
                      help="specify the term you want to use - GO, PFAM, PANTHER, KOG, KEGG, ALL")
    parser.add_option("-c", "--config", dest="configfile", default="annotation/goHyperG.config",
                      help="specify the configuration file")
    (options, args) = parser.parse_args()
    return options, parser

##################################################
# load association files
##################################################

def loadAssociation(associationfile):
    fh = open(associationfile, 'r')
    allterms = {}
    allgenes = {}

    for i in fh.readlines():
        linesplit = i.split()
        if len(linesplit) < 2:
            continue

        gene = linesplit[0]
        term = linesplit[1]
        eachterm = term.split(",")

        for e in eachterm:
            if gene in allgenes.keys():
                if e in allgenes[gene]:
                    pass
                else:
                    allgenes[gene].append(e)
            else:
                allgenes[gene]=[]
                allgenes[gene].append(e)

            if e in allterms.keys():
                if gene in allterms[e]:
                    pass
                else:
                    allterms[e].append(gene)
            else:
                allterms[e]=[]
                allterms[e].append(gene)

    return allgenes,allterms

##################################################
# load gene list
##################################################

def loadGeneList(genelistfile):
    fh = open(genelistfile, 'r')
    genelist = []

    for i in fh.readlines():
        linesplit = i.split()
        gene = linesplit[0]

        if gene in genelist:
            pass
        else:
            genelist.append(gene)

    return genelist

##################################################
# load association name
##################################################

def loadassocname(descriptionfile):
    fh = open(descriptionfile, 'r')
    goterm = {}

    for i in fh.readlines():
        linesplit = i.split()
        goterm[linesplit.pop(0)] = " ".join(linesplit)

    return goterm

##################################################
# load association name
##################################################

def loadConfig(configfile, species, term):
    fh = open(configfile, 'r')

    associationfile=[]
    descriptionfile=[]
    for i in fh.readlines():
        linesplit = i.strip("\n").split("|")
        if linesplit[0] == species:
            if linesplit[1] == term:
                associationfile.append(linesplit[2])
                descriptionfile.append(linesplit[3])
            elif term == "ALL":
                associationfile.append(linesplit[2])
                descriptionfile.append(linesplit[3])

    return associationfile, descriptionfile


##################################################
# do hyperg test
##################################################

def doHyperG(genelist, allgenes, allterms, assocname):

    geneswithterms = allgenes.keys()
    termswithgenes = allterms.keys()

    M=len(geneswithterms)
    N=len(list(set(geneswithterms).intersection(set(genelist))))

    pvalues=[]
    termsingenelist=[]
    termsinbackground=[]
    termname=[]

    for t in termswithgenes:
        n = len(allterms[t])
        x = len(list(set(allterms[t]).intersection(set(genelist))))
        if x == 0:
            continue

        pvalue = 1.0 - hypergeom.cdf(x,M,n,N)
        pvalues.append(pvalue)

        termsingenelist.append(x)
        termsinbackground.append(n)
        termname.append(t)

    adjpvalue = list(fdrcorrection0(pvalues)[1])

    print("\t".join(["Term annotation", "pvalue", "fdr adj pvalue","Background","Expected","GeneList","Observed","Genes"]))
    for u in range(0,len(adjpvalue)):
        gotermname = termname[u]
        if termname[u] in assocname.keys():
            gotermname = assocname[termname[u]]
        print("\t".join([gotermname,
                         str(pvalues[u]),
                         str(adjpvalue[u]),
                         str(M),
                         str(termsinbackground[u]),
                         str(N),
                         str(termsingenelist[u]),
                         ",".join(list(set(allterms[termname[u]]).intersection(set(genelist))))]
                        )
              )



##################################################
##### MAIN #######################################
##################################################

if __name__ == '__main__':
    allarguments,parser = parseArguments()

    # exit if no args provided
    if allarguments.species is None or allarguments.term is None:
        parser.print_help()
        sys.exit(1)

    associationfile, descriptionfile = loadConfig(allarguments.configfile, allarguments.species, allarguments.term)

    for i in range(len(associationfile)):

        allgenes, allterms = loadAssociation(str(associationfile[i]))
        genelist = loadGeneList(str(allarguments.genelistfile))
        assocname = loadassocname(str(descriptionfile[i]))
        doHyperG(genelist, allgenes, allterms, assocname)

