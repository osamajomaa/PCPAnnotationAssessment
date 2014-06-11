#!/usr/bin/env python
from collections import OrderedDict
from Bio import Entrez as ez
import unicodedata
import operator
import cPickle
import sys
import io
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
ez.email = "jomaao@miamioh.edu"

def uniquify_list(crosspec_rec):
    
    visited = []
    new_crosspec_rec = []
    for pair in crosspec_rec:
        if pair['pmid'] not in visited:
            new_crosspec_rec.append(pair)
            visited.append(pair['pmid'])
    return new_crosspec_rec

def mostcited_pmids(crosspec_pmid, spec1, spec2, top, paper_type, qualifier):
    
    if qualifier:
        s1s1_data, s1s2_data, s1n_data = get_mostcited_nonmethods_data(crosspec_pmid, spec1, spec2)
    else:
        s1s1_data, s1s2_data, s1n_data = get_mostcited_data(crosspec_pmid, spec1, spec2)
    
    build_mostcited_file(s1s1_data, spec1, spec1, top, paper_type, qualifier)
    build_mostcited_file(s1s2_data, spec1, spec2, top, paper_type, qualifier)
    build_mostcited_file(s1n_data, spec1, 'Others', top, paper_type, qualifier)


def build_mostcited_file(data, spec1, spec2, top, paper_type, qualifier):
    
    if len(data) < top:
        top = len(data)
    mostcited_pmid = OrderedDict(data.items()[:top])
    file_contents = "Top "+str(top)+" "+spec2+" papers that are cited by "+spec1+" papers:\n\n"
    for pmid in mostcited_pmid:
        handle = ez.efetch(db="pubmed", id=[pmid], rettype="medline", retmode="xml")
        records = ez.parse(handle)
        for record in records:
            file_contents += "PMID: " + pmid + "\n"
            file_contents += "Number of times this paper has been cited = " + str(mostcited_pmid[pmid]) + "\n"
            
            title = ""
            try:
                title = record['MedlineCitation']['Article']['ArticleTitle']
            except Exception:
                pass        
            file_contents += "Title: " + title +"\n"
            
            journal_title = ""
            try:
                journal_title = record['MedlineCitation']['Article']['Journal']['Title']
            except Exception:
                pass        
            file_contents += "Journal: " + journal_title +"\n"
            
            pubyear = ""
            try:
                pubyear = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
            except Exception:
                pass     
            pubmonth = ""
            try:
                pubmonth = record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Month']
            except Exception:
                pass
            file_contents += "Publication Date: " + pubmonth + ', ' + pubyear +"\n"
            
            abstract = ""
            try:
                for abst in record['MedlineCitation']['Article']['Abstract']['AbstractText']:
                    abstract += abst
            except Exception:
                pass        
            file_contents += "Abstract:\n" + abstract +"\n"        
            
            file_contents += "\n\n\n\n\n\n"
    if qualifier:
        file = open(os.path.join(CURR_PATH,"Ref_Catgzn/"+paper_type+"/"+spec1+"/nonmethods/top"+spec2+"_citedby_"+spec1+".txt"), 'w')
    else:
        file = open(os.path.join(CURR_PATH,"Ref_Catgzn/"+paper_type+"/"+spec1+"/all/top"+spec2+"_citedby_"+spec1+".txt"), 'w')
    #file_contents = unicodedata.normalize('NFKD', file_contents).encode('ascii','ignore')
    file.write(file_contents)
    file.close()


def get_mostcited_data(crosspec_pmid, spec1, spec2):
    s1s1_hist_data = OrderedDict()
    s1s2_hist_data = OrderedDict()
    s1n_hist_data = OrderedDict()
    for pmid in crosspec_pmid.keys():
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['species'] == spec1:
                if pair['pmid'] not in s1s1_hist_data:
                    s1s1_hist_data[pair['pmid']] = 1
                else:
                    s1s1_hist_data[pair['pmid']] += 1
            elif pair['species'] == spec2:
                if pair['pmid'] not in s1s2_hist_data:
                    s1s2_hist_data[pair['pmid']] = 1
                else:
                    s1s2_hist_data[pair['pmid']] += 1
            elif pair['species'] == 'None':
                if pair['pmid'] not in s1n_hist_data:
                    s1n_hist_data[pair['pmid']] = 1
                else:
                    s1n_hist_data[pair['pmid']] += 1
            elif pair['species'] == spec1+","+spec2:
                if pair['pmid'] not in s1s1_hist_data:
                    s1s1_hist_data[pair['pmid']] = 1
                else:
                    s1s1_hist_data[pair['pmid']] += 1
                if pair['pmid'] not in s1s2_hist_data:
                    s1s2_hist_data[pair['pmid']] = 1
                else:
                    s1s2_hist_data[pair['pmid']] += 1
                
    s1s1_hist_data = OrderedDict(sorted(s1s1_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    s1s2_hist_data = OrderedDict(sorted(s1s2_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    s1n_hist_data = OrderedDict(sorted(s1n_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    return s1s1_hist_data, s1s2_hist_data, s1n_hist_data 


def get_mostcited_nonmethods_data(crosspec_pmid, spec1, spec2):
    s1s1_hist_data = OrderedDict()
    s1s2_hist_data = OrderedDict()
    s1n_hist_data = OrderedDict()
    for pmid in crosspec_pmid.keys():
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['qualifier'] == "False":
                if pair['species'] == spec1:
                    if pair['pmid'] not in s1s1_hist_data:
                        s1s1_hist_data[pair['pmid']] = 1
                    else:
                        s1s1_hist_data[pair['pmid']] += 1
                elif pair['species'] == spec2:
                    if pair['pmid'] not in s1s2_hist_data:
                        s1s2_hist_data[pair['pmid']] = 1
                    else:
                        s1s2_hist_data[pair['pmid']] += 1
                elif pair['species'] == 'None':
                    if pair['pmid'] not in s1n_hist_data:
                        s1n_hist_data[pair['pmid']] = 1
                    else:
                        s1n_hist_data[pair['pmid']] += 1
                elif pair['species'] == spec1+","+spec2:
                    if pair['pmid'] not in s1s1_hist_data:
                        s1s1_hist_data[pair['pmid']] = 1
                    else:
                        s1s1_hist_data[pair['pmid']] += 1
                    if pair['pmid'] not in s1s2_hist_data:
                        s1s2_hist_data[pair['pmid']] = 1
                    else:
                        s1s2_hist_data[pair['pmid']] += 1
                
    s1s1_hist_data = OrderedDict(sorted(s1s1_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    s1s2_hist_data = OrderedDict(sorted(s1s2_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    s1n_hist_data = OrderedDict(sorted(s1n_hist_data.items(), key=operator.itemgetter(1), reverse=True))
    return s1s1_hist_data, s1s2_hist_data, s1n_hist_data 

def name_mapper(name):
    if name == 'human':
        return 'Humans'
    elif name == 'mouse':
        return 'Mice'

def initializer(spec1, spec2, paper_type, top, mflag):
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"Ref_Catgzn/"+paper_type+"/"+spec1+"_citations")))
    mostcited_pmids(crosspec_pmid, name_mapper(spec1), name_mapper(spec2), int(top), paper_type, mflag)
    
    
if __name__ == "__main__":

#     spec1 = sys.argv[1]
#     spec2 = sys.argv[2]
#     paper_type = sys.argv[3]
#     top = sys.argv[4]
#     mflag = True if sys.argv[5] == 'f' else False 
    initializer('human', 'mouse', 'nonreview', 50, False)
    initializer('human', 'mouse', 'nonreview', 50, True)
    initializer('human', 'mouse', 'review', 50, False)
    initializer('human', 'mouse', 'review', 50, True)
    
    initializer('mouse', 'human', 'nonreview', 50, False)
    initializer('mouse', 'human', 'nonreview', 50, True)
    initializer('mouse', 'human', 'review', 50, False)
    initializer('mouse', 'human', 'review', 50, True)
    
    