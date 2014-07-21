from Bio.UniProt import GOA
import httplib
from urllib import urlencode
import json
import os
import cPickle
from _mysql_exceptions import Error

RESPONSE_FORMAT = "application/json"
API_KEY = "eeeeac63bdcbd8551deacd2d4d445a00"
CONN_STRING = "api.elsevier.com:80"
CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def get_pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contained.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    unigoa_file = open(gaf_file)
    pmids = {}
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                pmids[pmid] = None
        
    return list(pmids.keys())

def scopusID_to_pmID(scopus_id):
    
    conn = httplib.HTTPConnection(CONN_STRING)
    params = urlencode({"field": "pubmed-id"})
    conn.putrequest("GET", "/content/abstract/scopus_id/"+scopus_id+"?"+params)
    conn.putheader("Accept", RESPONSE_FORMAT)    
    conn.putheader("X-ELS-APIKey", API_KEY)    
    conn.endheaders()
    response = conn.getresponse()
    print scopus_id
    refObj = json.loads(response.read())
    status = True
    payload = (response.status, response.reason)
    if response.status != 200:
        status = False
    else:
        try:
            payload = refObj["abstracts-retrieval-response"]["coredata"]["pubmed-id"]
        except Exception:
            status = False
            payload = (100, "Not a pubmed article")
    return status, payload

def get_ref_list(biblgy):
    
    query = ""
    for ref in biblgy:
        if 'scopus-id' in ref:
            query += "scopus-id%28" + ref['scopus-id'] + "%29" + "+OR+"
    
    query = query[:-4]
    conn = httplib.HTTPConnection(CONN_STRING)
    conn.putrequest("GET", "/content/search/scopus?query="+query)
    conn.putheader("Accept", RESPONSE_FORMAT)    
    conn.putheader("X-ELS-APIKey", API_KEY)
    conn.endheaders()
    response = conn.getresponse()
    data = json.loads(response.read())
    references = []
    try:
        pmid_data = data['search-results']['entry']
    except Exception:
        pass
    else:
        for entry in pmid_data:
            if "pubmed-id" in entry:
                references.append(entry["pubmed-id"])
    return references

def ref_grabber(pmid_list):
    
    pmid_refs = {}
    paper_errors = {}
    serial = 1
    conn = httplib.HTTPConnection(CONN_STRING)
    for pmid in pmid_list:
        pmid_refs[pmid] = []
        conn.putrequest("GET", "/content/abstract/pubmed_id/"+pmid+"?view=REF")
        conn.putheader("Accept", RESPONSE_FORMAT)    
        conn.putheader("X-ELS-APIKey", API_KEY)   
        conn.endheaders()
        response = conn.getresponse()
        data = json.loads(response.read())
        if response.status != 200:
            paper_errors[pmid] = (str(serial), str(response.status), str(response.reason))
        else:
            try:
                refs = data['abstracts-retrieval-response']['references']['reference']
            except Exception:
                print pmid
                paper_errors[pmid] = (str(serial), "100", "No references")
            else:
                pmid_refs[pmid] = get_ref_list(refs)
        serial += 1
                                             
    return pmid_refs, paper_errors

def create_ref_file(pmid_refs, species):
    fHandler = open("SCOPUS_API_Files/"+"pmid_pmid_"+species, 'wb')
    cPickle.dump(pmid_refs, fHandler)
    fHandler.close()

def create_err_file(paper_errors, species):
    err_file = "Number of papers not found in Scopus = " + str(len(paper_errors)) + "\n\n"
    for paper in paper_errors:
        err_file += paper_errors[paper][0] + "    "+ paper + ": (" + paper_errors[paper][1] + ", " + paper_errors[paper][2] + ")\n" 
    file = open("Erronous_Results/"+"not_found_scopus_pmid_"+species+".txt", 'wb')
    file.write(err_file)
    file.close()
    
if __name__ == "__main__":
    
    pmid_list = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/mouse/pmid_mouse")))
    pmid_refs, paper_errors = ref_grabber(pmid_list[:4])
    create_ref_file(pmid_refs, "mouse")
    create_err_file(paper_errors, "mouse")
    print "finished Mouse!"
    
#     pmid_list = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_human")))
#     pmid_refs, paper_errors = ref_grabber(pmid_list)
#     create_ref_file(pmid_refs, "human")
#     create_err_file(paper_errors, "human")
    
