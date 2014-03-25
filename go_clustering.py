from __future__ import division
from collections import OrderedDict
from operator import itemgetter
from sets import Set
from Bio.UniProt import GOA
import cPickle
import os

import sys

CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def build_clusters(species):
    """
        Build GO Clusters from a species gene association file. The cluster contains a representative GO term, 
        the proteins annotated to this term and all the papers that those proteins appear in.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    go_clusters = {}
    pmid_go = {}
    pmid_prot = {}
    unigoa_file = open(os.path.join(CURR_PATH,"GOA_Files/gene_association.goa_"+species))
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                if inrec['GO_ID'] not in go_clusters:
                    go_clusters[inrec['GO_ID']] = {'proteins':Set([inrec['DB_Object_ID']]), 'papers':Set([pmid])}
                else:
                    go_clusters[inrec['GO_ID']]['proteins'].add(inrec['DB_Object_ID'])
                    go_clusters[inrec['GO_ID']]['papers'].add(pmid)
                if pmid not in pmid_go:
                    pmid_go[pmid] = Set([inrec['GO_ID']])
                else:
                    pmid_go[pmid].add(inrec['GO_ID'])
                if pmid not in pmid_prot:
                    pmid_prot[pmid] = Set([inrec['DB_Object_ID']])
                else:
                    pmid_prot[pmid].add(inrec['DB_Object_ID'])
    
    pickle_data(go_clusters, os.path.join(CURR_PATH, "Pickled_Data/go_clusters_"+species))
    pickle_data(pmid_go, os.path.join(CURR_PATH, "Pickled_Data/pmid_go_"+species))
    pickle_data(pmid_prot, os.path.join(CURR_PATH, "Pickled_Data/pmid_prot_"+species))


def pickle_data(data, path):
    
    fHandler = open(path, 'wb')
    cPickle.dump(data, fHandler)
    fHandler.close()

def load_data(species):
    go_clusters = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/go_clusters_"+species)))
    pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/pmid_go_"+species)))
    pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/pmid_prot_"+species)))
    
    return go_clusters, pmid_go, pmid_prot

    
def get_paper_cohesion(go_clusters):
    
    goterm_per_paper = {}
    paper_cohesion_ratio = OrderedDict()
    for go_term in go_clusters.keys():
        for paper in go_clusters[go_term]['papers']:
            if paper not in goterm_per_paper:
                goterm_per_paper[paper] = 1
            else:
                goterm_per_paper[paper] += 1
    
    total_terms = len(go_clusters)
    for paper in goterm_per_paper.keys():
        paper_cohesion_ratio[paper] = ((total_terms-goterm_per_paper[paper])*100)/total_terms
    x = OrderedDict(sorted(paper_cohesion_ratio.items(), key=itemgetter(1), reverse=True))
    
    return x

def top_cohesive_papers(top, paper_cohesion_ratio, pmid_prot):
    
    
    #top_pmids = OrderedDict(paper_cohesion_ratio.items()[:top])
    i = 0
    top_pmids = OrderedDict()
    for pmid in paper_cohesion_ratio:
        if i >= top:
            break
        if len(pmid_prot[pmid]) >= 4:
            top_pmids[pmid] = paper_cohesion_ratio[pmid]
            i += 1
            
    output = "{0[0]:<15}{0[1]:<15}".format("PMID Cohesion_Rate".split()) + "\n"
    line = [0,0]
    for pmid in top_pmids.keys():
        line[0] = pmid
        line[1] = str(top_pmids[pmid])
        output += "{0[0]:<15}{0[1]:<15}".format(line) + "\n"
    
    return output

def top_cohesive_papers_protsAndterms(top, paper_cohesion_ratio, pmid_prot, pmid_go):
    
    #top_pmids = OrderedDict(paper_cohesion_ratio.items()[:top])
    i = 0
    top_pmids = OrderedDict()
    for pmid in paper_cohesion_ratio:
        if i >= top:
            break
        if len(pmid_prot[pmid]) >= 4:
            top_pmids[pmid] = paper_cohesion_ratio[pmid]
            i += 1
            
    output = ""
    for pmid in top_pmids.keys():
        output += "PMID = " + pmid + "\n"
        output += "Cohesion Rate = " + str(top_pmids[pmid]) + "\n"
        output += "Proteins:\n"
        for protein in pmid_prot[pmid]:
            output += protein + ", "
        output = output[:-2]
        output += "\n"
        output += "GO Terms:\n"
        for go_term in pmid_go[pmid]:
            output += go_term + ", "
        output = output[:-2]
        output += "\n\n"
    
    return output


def print_output(output):
    sys.stdout.write(output)


if __name__ == "__main__":

    #build_clusters("dicty")
    go_clusters, pmid_go, pmid_prot = load_data("dicty")
    paper_cohesion_ratio = get_paper_cohesion(go_clusters)
    output = top_cohesive_papers_protsAndterms(20, paper_cohesion_ratio, pmid_prot, pmid_go)
    #output = top_cohesive_papers(10, paper_cohesion_ratio)
    print output

    
