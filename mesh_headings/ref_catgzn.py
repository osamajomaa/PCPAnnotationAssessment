#!/usr/bin/env python
from Bio import Entrez as ez
import cPickle
import json
import sys
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
ez.email = "jomaao@miamioh.edu"

def search_species(meshHeadings):
    desc_names = []
    for heading in meshHeadings:
        for qual_name in heading['QualifierName']:
            desc_names.append(qual_name)
        desc_names.append(heading['DescriptorName'])
    return desc_names

def search_headings(mesh_hds, term):
    for heading in mesh_hds:
        if heading.lower().find(term) != -1:
            return True
    return False

def add_to_crosspec(crosspec_pmid, pubmed_val, desc_names, pmid):
    if 'Humans' in desc_names and 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans,Mice'})
    elif 'Humans' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans'})
    elif 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Mice'})
    else:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'None'})
    return crosspec_pmid


def pickle_data(data, path):
    
    fHandler = open(path, 'wb')
    cPickle.dump(data, fHandler)
    fHandler.close()      

def medline(pmid_pmid, spec1, spec2):
    
    #crosspec_pmid = {}
    crosspec_pmid_methods = {}
    crosspec_pmid_review = {}
    crosspec_pmid_research = {}
    
    for pmid in pmid_pmid.keys():
        try:
            handle = ez.efetch(db="pubmed", id=pmid_pmid[pmid], rettype="medline", retmode="xml")
            crosspec_pmid_methods[pmid] = []
            crosspec_pmid_review[pmid] = []
            crosspec_pmid_research[pmid] = []
            records = ez.parse(handle)
            
            for pubmed_rec in records:
                desc_names = search_species(pubmed_rec['MedlineCitation']['MeshHeadingList'])
                if search_headings(desc_names, 'methods'):
                    crosspec_pmid_methods = add_to_crosspec(crosspec_pmid_methods, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)
                elif search_headings(desc_names, 'review'):
                    crosspec_pmid_review = add_to_crosspec(crosspec_pmid_review, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)
                else:
                    crosspec_pmid_research = add_to_crosspec(crosspec_pmid_research, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)                    
        except Exception:
            pass
    
    pickle_data(crosspec_pmid_methods, os.path.join(CURR_PATH, "mesh_headings/methods/"+spec1+"_citations"))
    pickle_data(crosspec_pmid_review, os.path.join(CURR_PATH, "mesh_headings/review/"+spec1+"_citations"))
    pickle_data(crosspec_pmid_research, os.path.join(CURR_PATH, "mesh_headings/research/"+spec1+"_citations"))

if __name__ == "__main__":
    
    main_spec = sys.argv[1]
    sec_spec = sys.argv[2]
    #main_spec = 'mouse'
    #sec_spec = 'human'
    pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+"pmid_pmid_"+main_spec)))
    medline(pmid_pmid, main_spec, sec_spec)
    print "Done with" + main_spec
