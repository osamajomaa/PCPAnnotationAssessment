#!/usr/bin/env python
import cPickle
import sys
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def uniquify_list(crosspec_rec):
    
    visited = []
    new_crosspec_rec = []
    for pair in crosspec_rec:
        if pair['pmid'] not in visited:
            new_crosspec_rec.append(pair)
            visited.append(pair['pmid'])
    return new_crosspec_rec

def build_crosspec_files(crosspec_pmid, spec1, spec2, paper_type):
    desc_text_s1s1 = "{0[0]:<15}{0[1]:<15}".format([spec1, spec1])+ "\n"
    desc_text_s1s2 = "{0[0]:<15}{0[1]:<15}".format([spec1, spec2])+ "\n"
    desc_text_s1n = "{0[0]:<15}{0[1]:<15}".format([spec1, 'Other'])+ "\n"
    count_s1s1 = 0
    count_s1s2 = 0
    count_s1n = 0
    for pmid in crosspec_pmid.keys():
        cit_s1s1 = ""
        cit_s1s2 = ""
        cit_s1n = ""
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['species'] == spec1:
                cit_s1s1 += pair['pmid'] + ", "
                count_s1s1 += 1
            elif pair['species'] == spec2:
                cit_s1s2 += pair['pmid'] + ", "
                count_s1s2 += 1
            elif pair['species'] == 'None':
                cit_s1n += pair['pmid'] + ", "
                count_s1n += 1
            elif pair['species'] == spec1+","+spec2:
                cit_s1s2 += pair['pmid'] + ", " 
                cit_s1s1 += pair['pmid'] + ", "
                count_s1s1 += 1
                count_s1s2 += 1        
        desc_text_s1s1 += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1s1[:-2]]) + "\n"
        desc_text_s1s2 += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1s2[:-2]]) + "\n"
        desc_text_s1n += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1n[:-2]]) + "\n"
    desc_text_s1s1 = "Number of " + spec1 + " to " + spec1 + " citations = " + str(count_s1s1) + "\n\n" + desc_text_s1s1
    desc_text_s1s2 = "Number of " + spec1 + " to " + spec2 + " citations = " + str(count_s1s2) + "\n\n" + desc_text_s1s2
    desc_text_s1n  = "Number of " + spec1 + " to " + "Other" + " citations = " + str(count_s1n) + "\n\n" + desc_text_s1n
    
    desc_file_s1s1 = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_"+spec1), 'w')
    desc_file_s1s1.write(desc_text_s1s1)
    desc_file_s1s2 = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_"+spec2), 'w')
    desc_file_s1s2.write(desc_text_s1s2)
    desc_file_s1n = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_Other"), 'w')
    desc_file_s1n.write(desc_text_s1n)

def name_mapper(name):
    if name == 'human':
        return 'Humans'
    elif name == 'mouse':
        return 'Mice'
    
if __name__ == "__main__":

    spec1 = sys.argv[1]
    spec2 = sys.argv[2]
    paper_type = sys.argv[3]
    
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"_citations")))
    build_crosspec_files(crosspec_pmid, name_mapper(spec1), name_mapper(spec2), paper_type)
    
    