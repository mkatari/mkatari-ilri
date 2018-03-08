"""
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
"""
__author__ = 'manpreetkatari'

import os.path as path
from argparse import ArgumentParser
from collections import defaultdict

from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection

BASE_DIR = path.dirname(path.realpath(__file__))


##################################################
# Do all the parsing of command line options here
##################################################

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("genelistfile", default="cassava_genelist.txt",
                        help="specify path to file with list of genes")
    parser.add_argument("-s", "--species", dest="species", default="cassavaV6", choices=("cassavaV6",),
                        help="specify the species to use")
    parser.add_argument("-t", "--term", dest="term", default="ALL",
                        choices=("GO", "PFAM", "PANTHER", "KOG", "KEGG", "ALL"),
                        help="specify the term you want to use - GO, PFAM, PANTHER, KOG, KEGG, ALL")
    parser.add_argument("-c", "--config", dest="configfile",
                        default=path.join(BASE_DIR, "annotation/goHyperG.config"),
                        help="specify the configuration file")
    return parser.parse_args()


##################################################
# load association files
##################################################

def load_association(association_file):
    all_terms = defaultdict(set)
    all_genes = defaultdict(set)

    with open(association_file, 'r') as fh:
        for line in fh:
            try:
                gene, term = line.split()

                for e in term.split(","):
                    all_genes[gene].add(e)
                    all_terms[e].add(gene)
            except ValueError:
                pass

    return all_genes, all_terms


##################################################
# load gene list
##################################################

def load_gene_list(genelistfile):
    with open(genelistfile, 'r') as fh:
        return {i.split()[0] for i in fh}


##################################################
# load association name
##################################################

def load_assoc_name(description_file):
    with open(description_file, 'r') as fh:
        return {name: desc for name, sep, desc in (line.rstrip("\n").partition("\t") for line in fh)}


##################################################
# load association name
##################################################

def load_config(config_file, species, term):
    with open(config_file, 'r') as fh:
        association_files = []
        description_files = []
        for line in fh:
            s, t, af, df = line.rstrip("\n").split("|")
            if s == species:
                if t == term:
                    association_files.append(path.join(BASE_DIR, af))
                    description_files.append(path.join(BASE_DIR, df))
                elif term == "ALL":
                    association_files.append(path.join(BASE_DIR, af))
                    description_files.append(path.join(BASE_DIR, df))

        return association_files, description_files


##################################################
# do hyperg test
##################################################

def do_hyper_geom(genelist, allgenes, allterms, assocname):
    M = len(allgenes)
    N = len(all_genes.keys() & genelist)

    pvalues = []
    termsingenelist = []
    termsinbackground = []
    termname = []

    for t, genes in allterms.items():
        x = len(genes & genelist)
        if not x:
            continue
        n = len(genes)

        pvalues.append(1.0 - hypergeom.cdf(x, M, n, N))

        termsingenelist.append(x)
        termsinbackground.append(n)
        termname.append(t)

    adjpvalue = fdrcorrection(pvalues)[1]

    print("\t".join(
        ["Term annotation", "pvalue", "fdr adj pvalue", "Background", "Expected", "GeneList", "Observed", "Genes"]))
    for p, adj_p, tb, tl, tn in zip(pvalues, adjpvalue, termsinbackground, termsingenelist, termname):
        try:
            gotermname = tn + " " + assocname[tn]
        except KeyError:
            gotermname = tn

        print("\t".join([gotermname,
                         str(p),
                         str(adj_p),
                         str(M),
                         str(tb),
                         str(N),
                         str(tl),
                         ",".join(allterms[tn] & genelist)]
                        )
              )


##################################################
#  MAIN
##################################################

if __name__ == '__main__':

    args = parse_arguments()

    association_files, description_files = load_config(args.configfile, args.species, args.term)
    genelist = load_gene_list(args.genelistfile)

    for af, df in zip(association_files, description_files):
        all_genes, all_terms = load_association(af)
        assocname = load_assoc_name(df)
        do_hyper_geom(genelist, all_genes, all_terms, assocname)
