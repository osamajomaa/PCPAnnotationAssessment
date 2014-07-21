#!/usr/bin/env python
from collections import OrderedDict
import operator
import cPickle
import pylab
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

def draw_cited_histograms(crosspec_pmid, spec1, spec2, top, paper_type):
    
    s1s1_data, s1s2_data, s1n_data = get_cited_histogram_data(crosspec_pmid, spec1, spec2)
    
    draw_cited_histogram_helper(s1s1_data, spec1, spec1, top, 'red', paper_type)
    draw_cited_histogram_helper(s1s2_data, spec1, spec2, top, 'green', paper_type)
    draw_cited_histogram_helper(s1n_data, spec1, 'Others', top, 'blue', paper_type)
    

def draw_cited_histogram_helper(hist_data, spec1, spec2, top, bcolor, paper_type):
    #print hist_data.values()
    if len(hist_data) < top:
        top = len(hist_data)
    xvalues = hist_data.keys()[:top]
    yaxis = hist_data.values()[:top]
    xaxis = list(range(0,top))
    
    pylab.bar(xaxis, yaxis, align='center', color=bcolor)
    pylab.gcf().subplots_adjust(bottom=0.25)    
    pylab.xticks(xaxis, xvalues, rotation=90)
    

    pylab.title("Top "+str(top)+" "+spec2+" papers that are cited by "+spec1+" papers")

    pylab.xlabel(spec2 + " Papers")
    pylab.ylabel("Number of times the paper has been cited by " + spec1 + " papers")        
    pylab.grid(True)
    pylab.tight_layout()
    pylab.savefig(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/histograms/top"+spec2+"_citedby_"+spec1+".png"))
    pylab.close()


def get_cited_histogram_data(crosspec_pmid, spec1, spec2):
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

def name_mapper(name):
    if name == 'human':
        return 'Humans'
    elif name == 'mouse':
        return 'Mice'
    
if __name__ == "__main__":

    spec1 = sys.argv[1]
    spec2 = sys.argv[2]
    paper_type = sys.argv[3]
    top = sys.argv[4]
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"_citations")))
    draw_cited_histograms(crosspec_pmid, name_mapper(spec1), name_mapper(spec2), int(top), paper_type)
    
    