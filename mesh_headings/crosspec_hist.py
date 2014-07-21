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

def draw_histograms(crosspec_pmid, spec1, spec2, paper_type):
    
    s1s1_data, s1s2_data, s1n_data = get_histogram_data(crosspec_pmid, spec1, spec2)
    
    draw_histogram_helper(s1s1_data, spec1, spec1, 'red', paper_type)
    draw_histogram_helper(s1s2_data, spec1, spec2, 'green', paper_type)
    draw_histogram_helper(s1n_data, spec1, 'Others', 'blue', paper_type)
    
    
def draw_histogram_helper(hist_data, spec1, spec2, bcolor, paper_type):
    
    xaxis = []
    yaxis = []
    for key,value in hist_data.items():
        xaxis.append(value[0])
        yaxis.append(key)
    
    pylab.bar(xaxis, yaxis, align='center', color=bcolor)
    pylab.gcf().subplots_adjust(bottom=0.25)    
    pylab.xticks(xaxis, xaxis, rotation=90)
    

    pylab.title("Distribution of "+spec1+"-"+spec2 + " Citations")

    pylab.xlabel("Number of citations in " + spec2)
    pylab.ylabel("Number of papers in " + spec1)        
    pylab.grid(True)
    pylab.tight_layout()
    pylab.savefig(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/histograms/"+spec1+"_cit_"+spec2+".png"))
    pylab.close()
    
    
def get_histogram_data(crosspec_pmid, spec1, spec2):
    
    s1s1_hist_data = {}
    s1s2_hist_data = {}
    s1n_hist_data = {}
    s1s1_total = 0
    s1s2_total = 0
    s1n_total = 0
    for pmid in crosspec_pmid.keys():
        count_s1s1 = 0
        count_s1s2 = 0
        count_s1n = 0
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['species'] == spec1:
                count_s1s1 += 1
            elif pair['species'] == spec2:
                count_s1s2 += 1
            elif pair['species'] == 'None':
                count_s1n += 1
            elif pair['species'] == spec1+","+spec2:
                count_s1s2 += 1
                count_s1s1 += 1
        if count_s1s1 > 0:
            if count_s1s1 not in s1s1_hist_data:
                s1s1_hist_data[count_s1s1] = [1,[]]
            else:
                s1s1_hist_data[count_s1s1][0] += 1
                s1s1_hist_data[count_s1s1][1].append(pair['pmid'])
        if count_s1s2 > 0:
            if count_s1s2 not in s1s2_hist_data:
                s1s2_hist_data[count_s1s2] = [1,[]]
            else:
                s1s2_hist_data[count_s1s2][0] += 1
                s1s2_hist_data[count_s1s2][1].append(pair['pmid'])
        if count_s1n > 0:
            if count_s1n not in s1n_hist_data:
                s1n_hist_data[count_s1n] = [1,[]]
            else:
                s1n_hist_data[count_s1n][0] += 1
                s1n_hist_data[count_s1n][1].append(pair['pmid'])
        s1s1_total += count_s1s1
        s1s2_total += count_s1s2
        s1n_total += count_s1n
    
#     for c in s1s1_hist_data:
#         print s1s1_hist_data[c][0]
#         s1s1_hist_data[c][0] = round(s1s1_hist_data[c][0]*100.0/s1s1_total, 3)
#         print s1s1_hist_data[c][0]
#     for c in s1s2_hist_data:
#         s1s2_hist_data[c][0] = round(s1s2_hist_data[c][0]*100.0/s1s2_total, 3)
#     for c in s1n_hist_data:        
#         s1n_hist_data[c][0] = round(s1n_hist_data[c][0]*100.0/s1n_total, 3)    
                
    s1s1_hist_data = dict(sorted(s1s1_hist_data.items()))
    s1s2_hist_data = dict(sorted(s1s2_hist_data.items()))
    s1n_hist_data = dict(sorted(s1n_hist_data.items()))
    return s1s1_hist_data, s1s2_hist_data, s1n_hist_data 

                
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

    #spec1 = sys.argv[1]
    #spec2 = sys.argv[2]
    #paper_type = sys.argv[3]
    spec1 = 'human'
    spec2 = 'mouse'
    paper_type = 'methods'    
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"_citations")))
    draw_histograms(crosspec_pmid, name_mapper(spec1), name_mapper(spec2), paper_type)
    
    
