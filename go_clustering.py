from __future__ import division
from collections import OrderedDict
from operator import itemgetter
from sets import Set
from Bio.UniProt import GOA
import cPickle
import os
import pylab


CURR_PATH = os.path.dirname(os.path.realpath(__file__))
GO_Species = {'P':'bio_proc', 'F':'mol_func', 'C':'cell_comp'}

def add_to_pmid(pmid_x, pmid, x):
    
    if pmid not in pmid_x:
        pmid_x[pmid] = Set([x])
    else:
        pmid_x[pmid].add(x)

def build_clusters(species):
    """
        Build GO Clusters from a species gene association file. The cluster contains a representative GO term, 
        the proteins annotated to this term and all the papers that those proteins appear in.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    
    pmid_go_mf = OrderedDict()
    pmid_go_bp = OrderedDict()
    pmid_go_cc = OrderedDict()
    pmid_go = OrderedDict()
    
    pmid_prot_mf = OrderedDict()
    pmid_prot_bp = OrderedDict()
    pmid_prot_cc = OrderedDict()
    pmid_prot = OrderedDict()
    
    go_prot_mf = OrderedDict()
    go_prot_bp = OrderedDict()
    go_prot_cc = OrderedDict()
    go_prot = OrderedDict()
    
    
    unigoa_file = open(os.path.join(CURR_PATH,"GOA_Files/gene_association.goa_"+species))
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                if inrec['Aspect'] == 'P':
                    add_to_pmid(pmid_go_bp, pmid, inrec['GO_ID'])
                    add_to_pmid(pmid_prot_bp, pmid, inrec['DB_Object_ID'])
                    add_to_pmid(go_prot_bp, inrec['GO_ID'], inrec['DB_Object_ID'])
                elif inrec['Aspect'] == 'F':
                    add_to_pmid(pmid_go_mf, pmid, inrec['GO_ID'])
                    add_to_pmid(pmid_prot_mf, pmid, inrec['DB_Object_ID'])
                    add_to_pmid(go_prot_mf, inrec['GO_ID'], inrec['DB_Object_ID'])
                elif inrec['Aspect'] == 'C':
                    add_to_pmid(pmid_go_cc, pmid, inrec['GO_ID'])
                    add_to_pmid(pmid_prot_cc, pmid, inrec['DB_Object_ID'])
                    add_to_pmid(go_prot_cc, inrec['GO_ID'], inrec['DB_Object_ID'])
                add_to_pmid(pmid_go, pmid, inrec['GO_ID'])
                add_to_pmid(pmid_prot, pmid, inrec['DB_Object_ID'])
                add_to_pmid(go_prot, inrec['GO_ID'], inrec['DB_Object_ID'])
                
        
    pickle_data(pmid_go_mf, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_go_mf_"+species))
    pickle_data(pmid_go_cc, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_go_cc_"+species))
    pickle_data(pmid_go_bp, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_go_bp_"+species))
    pickle_data(pmid_go, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_go_"+species))
    
    pickle_data(pmid_prot_mf, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_prot_mf_"+species))
    pickle_data(pmid_prot_cc, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_prot_cc_"+species))
    pickle_data(pmid_prot_bp, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_prot_bp_"+species))
    pickle_data(pmid_prot, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/pmid_prot_"+species))
    
    pickle_data(go_prot_mf, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/go_prot_mf_"+species))
    pickle_data(go_prot_cc, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/go_prot_cc_"+species))
    pickle_data(go_prot_bp, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/go_prot_bp_"+species))
    pickle_data(go_prot, os.path.join(CURR_PATH, "Pickled_Data/"+species+"/go_prot_"+species))


def pickle_data(data, path):
    
    fHandler = open(path, 'wb')
    cPickle.dump(data, fHandler)
    fHandler.close()

def load_data(species, onto='all'):
    
    pmid_go = OrderedDict()
    pmid_prot = OrderedDict()
    
    if onto.lower() == 'all':
        pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_go_"+species)))
        pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_prot_"+species)))
        go_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/go_prot_"+species)))
        
    elif onto.lower() == 'mf':
        pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_go_mf_"+species)))
        pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_prot_mf_"+species)))
        go_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/go_prot_mf_"+species)))
    
    elif onto.lower() == 'cc':
        pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_go_cc_"+species)))
        pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_prot_cc_"+species)))
        go_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/go_prot_cc_"+species)))
    
    elif onto.lower() == 'bp':
        pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_go_bp_"+species)))
        pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_prot_bp_"+species)))
        go_prot = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/go_prot_bp_"+species)))
    
    return pmid_go, pmid_prot, go_prot

    

def calc_cohesion(pmid_go, pmid_prot, step, threshold=-1):
    
    '''
        First off, I group the papers (pmids) according to the number of GO terms they have.
        Ex: 1: pmid1, pmid; 2: pmid3, pmid4, pmid; ...
    '''
    go_cats = OrderedDict()
    
    for pmid in pmid_go.keys():
        tcount = len(pmid_go[pmid]) 
        if tcount not in go_cats:
            go_cats[tcount] = [pmid]
        else:
            go_cats[tcount].append(pmid)
    
    
    '''    
        Then, I start to categorize the papers according to the average number of GO terms per step
    '''
    go_cats =  OrderedDict(sorted(go_cats.items()))
    minm = 0
    maxm = step
    max_go = len(go_cats)
    avg_go_cat = []
    while maxm <= max_go:
        slc = OrderedDict(go_cats.items()[minm:maxm])
        summ = 0
        for go_num in slc.keys():
            summ += go_num
        avg = summ/step 
        pmids = []
        for pmid in slc.keys():
            pmids.extend(slc[pmid]) 
        print len(pmids)
        avg_go_cat.append((pmids, avg))
        minm = maxm
        maxm = maxm + step
    
    if minm < max_go:
        slc = OrderedDict(go_cats.items()[minm:maxm])
        summ = 0
        for go_num in slc.keys():
            summ += go_num
        avg = summ/len(slc)
        pmids = []
        for pmid in slc.keys():
            pmids.extend(slc[pmid]) 
        avg_go_cat.append((pmids, avg))
    
    
    '''
        Finally, I sort the papers (pmids) according to the number of GO terms and then
        according to the number of proteins removing the papers that have less than
        minimum number of proteins.
    '''
    print avg_go_cat
    pmid_compr = []
    for i in range(0, len(avg_go_cat)):
        go_range = OrderedDict()
        for j in range(0, len(avg_go_cat[i][0])):
            pmid = avg_go_cat[i][0][j]
            if len(pmid_prot[pmid]) > threshold:
                go_range[pmid] = len(pmid_prot[pmid])        
        go_range = OrderedDict(sorted(go_range.items(), key=itemgetter(1), reverse=True))
        print go_range
        pmid_compr.extend(list(go_range.keys()))
    
    return avg_go_cat, pmid_compr

def draw_go_pmid_histogram(go_cats, species, onto='all'):
        
    xaxis = []
    yaxis = []
    for (pmid_list, avg) in go_cats:
        xaxis.append(avg)
        yaxis.append(len(pmid_list))
        
    pylab.bar(xaxis, yaxis, align='center')
    pylab.gcf().subplots_adjust(bottom=0.25)    
    pylab.xticks(xaxis, xaxis, rotation=90)
    
    if onto.lower() == 'all':
        pylab.title("Number of papers per GO Category across all GO subontologies")
    elif onto.lower() == 'mf':
        pylab.title("Number of papers per GO Category across GO Molecular Function")      
    elif onto.lower() == 'bp':
        pylab.title("Number of papers per GO Category across GO Biological Process")
    elif onto.lower() == 'cc':
        pylab.title("Number of papers per GO Category across GO Cellular Component")

    pylab.xlabel("Average number of GO terms")
    pylab.ylabel('Number of PMIDs')        
    pylab.grid(True)
    pylab.tight_layout()
    pylab.savefig(os.path.join(CURR_PATH,"Cohesion/"+species+"/goPmid_"+species+"_"+onto+".png"))
    pylab.close()


if __name__ == "__main__":

    species = "dicty"
    onto = 'all'
    #build_clusters(species)
    pmid_go, pmid_prot, go_prot = load_data(species, onto=onto)
    go_cats, pmid_coh  = calc_cohesion(pmid_go, pmid_prot, 1)
    print pmid_coh
    draw_go_pmid_histogram(go_cats, species, onto)
    
    
    #go_clusters, pmid_go, pmid_prot = load_data("dicty")
    #paper_cohesion_ratio = get_paper_cohesion(go_clusters)
    #output = top_cohesive_papers_protsAndterms(20, paper_cohesion_ratio, pmid_prot, pmid_go)
    #output = top_cohesive_papers(10, paper_cohesion_ratio)
    #print output

    
