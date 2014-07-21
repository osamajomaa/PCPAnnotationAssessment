'''
    This script is to break down a list of requested pmids in a number of request sublists. 
    each sublist contains pmids that are to be used in a string to query SCOPUS API and get all
    the citations for these papers
'''

import cPickle
import sys
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def chunk_data(species, dest, numof_chunks):
    dum_pmids = cPickle.load(open(os.path.join(CURR_PATH,src)))
    set_size = len(dum_pmids)
    chunk_size = len(dum_pmids)/numof_chunks
    startof_chunk = 0
    endof_chunk = chunk_size if chunk_size <= set_size else set_size
    serial = 1
    if isinstance(dum_pmids, dict):
        while endof_chunk <= set_size:        
            chunk = dict(dum_pmids.items()[startof_chunk:endof_chunk])
            fhandler = open(os.path.join(CURR_PATH,dest+"/req"+str(serial)), 'wb')
            cPickle.dump(chunk, fhandler)
            fhandler.close()
            startof_chunk = endof_chunk
            endof_chunk = startof_chunk + chunk_size
            if 0 < (set_size-endof_chunk) < chunk_size:
                endof_chunk = set_size
            serial += 1
    else:        
        while endof_chunk <= set_size:
            
            chunk = dum_pmids[startof_chunk:endof_chunk]
            fhandler = open(os.path.join(CURR_PATH,dest+"/req"+str(serial)), 'wb')
            cPickle.dump(chunk, fhandler)
            fhandler.close()
            startof_chunk = endof_chunk
            endof_chunk = startof_chunk + chunk_size
            if 0 < (set_size-endof_chunk) < chunk_size:
                endof_chunk = set_size
            serial += 1



if __name__ == "__main__":
    
    #src = "Pickled_Data/"+species+"/pmid_"+species
    #dest = "SCOPUS_API_Files/"+species+"/Requests
    src = sys.argv[1]
    dest = sys.argv[2]
    numof_chunks = sys.argv[3]
    chunk_data(src, dest, int(numof_chunks))
    
    
    
    