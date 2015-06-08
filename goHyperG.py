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


##################################################
# Do all the parsing of command line options here
##################################################

def parseArguments():
    parser = OptionParser()
    parser.add_option("-g", "--genelist", dest="genelistfile", default="cassava_genelist.txt",
                      help="specify path to file with list of genes")
    parser.add_option("-a", "--associationfile", dest="associationfile", default="cassavaV61.go",
                      help="association file")
    parser.add_option("-d", "--description", dest="descriptionfile", default="gonames.txt",
                      help="file containing the description of the terms that are in association file")
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
 #       print(gene)
        term = linesplit[1]
  #      print(term)
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
 #       print(gene)

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
        print("\t".join([assocname[termname[u]],
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
    allarguments = parseArguments()

    # exit if no args provided
    if allarguments.genelistfile is None and allarguments.associationfile is None and allarguments.descriptionfile is None:
        parser.print_help()
        sys.exit(1)


    allgenes, allterms = loadAssociation(str(allarguments.associationfile))
    genelist = loadGeneList(str(allarguments.genelistfile))
    assocname = loadassocname(str(allarguments.descriptionfile))
    doHyperG(genelist, allgenes, allterms, assocname)

