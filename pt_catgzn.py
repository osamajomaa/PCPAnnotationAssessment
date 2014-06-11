#!/usr/bin/env python
from Bio import Entrez as ez
import cPickle
import json
import sys
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
ez.email = "jomaao@miamioh.edu"



def search_types(ptype_list, term):
    for ptype in ptype_list:
        if ptype.lower().find(term) != -1:
            return True
    return False

def search_species(meshHeadings):
    desc_names = []
    for heading in meshHeadings:
        for qual_name in heading['QualifierName']:
            desc_names.append(qual_name)
        desc_names.append(heading['DescriptorName'])
    return desc_names

def add_to_crosspec(crosspec_pmid, pubmed_val, desc_names, pmid, qualifier):
    if 'Humans' in desc_names and 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans,Mice', 'qualifier':str(qualifier)})
    elif 'Humans' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans', 'qualifier':str(qualifier)})
    elif 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Mice', 'qualifier':str(qualifier)})
    else:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'None', 'qualifier':str(qualifier)})
    return crosspec_pmid


def pickle_data(data, path):
    
    fHandler = open(path, 'wb')
    cPickle.dump(data, fHandler)
    fHandler.close()      

def medline(pmid_pmid, spec1, spec2):

    crosspec_pmid_review = {}
    crosspec_pmid_nonreview = {}
    
    for pmid in pmid_pmid.keys():
        try:
            handle = ez.efetch(db="pubmed", id=['9254694'], rettype="medline", retmode="xml")
            crosspec_pmid_review[pmid] = []
            crosspec_pmid_nonreview[pmid] = []
            records = ez.parse(handle)
            
            for pubmed_rec in records:
                types_list = pubmed_rec['MedlineCitation']['Article']['PublicationTypeList']
                desc_names = search_species(pubmed_rec['MedlineCitation']['MeshHeadingList'])
                qualifier = search_types(desc_names, "methods") 
                if search_types(types_list, 'review'):
                    crosspec_pmid_review = add_to_crosspec(crosspec_pmid_review, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid, qualifier)
                else:
                    crosspec_pmid_nonreview = add_to_crosspec(crosspec_pmid_nonreview, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid, qualifier)                    
        except Exception:
            pass
    
    pickle_data(crosspec_pmid_review, os.path.join(CURR_PATH, "Ref_Catgzn/review/"+spec1+"_citations"))
    pickle_data(crosspec_pmid_nonreview, os.path.join(CURR_PATH, "Ref_Catgzn/nonreview/"+spec1+"_citations"))

if __name__ == "__main__":
    
#     main_spec = sys.argv[1]
#     sec_spec = sys.argv[2]
    main_spec = 'mouse'
    sec_spec = 'human'
    pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+"pmid_pmid_"+main_spec)))
    medline(dict(pmid_pmid.items()), main_spec, sec_spec)
    print "Done with" + main_spec
