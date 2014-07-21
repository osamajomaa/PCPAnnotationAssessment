from Bio.UniProt import GOA
import httplib
from urllib import urlencode
import cPickle
import xml.etree.ElementTree as ET

RESPONSE_FORMAT = "text/xml"
API_KEY = "eeeeac63bdcbd8551deacd2d4d445a00"
CONN_STRING = "api.elsevier.com:80"

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
    root = ET.Element(ET.fromstring(response.read()))
    if response.status != 200:
        return False, (response.status, response.reason)
    else:
        for pubmed_id in root.tag.iter("{http://www.elsevier.com/xml/svapi/abstract/dtd}pubmed-id"):
            return True, pubmed_id.text
        return False, (100, "Not a pubmed article")
    
def get_ref_list(pmid_list):
    
    pmid_refs = {}
    ref_errors = {}
    paper_errors = {}
    conn = httplib.HTTPConnection(CONN_STRING)
    print len(pmid_list)
    i = 0
    for pmid in pmid_list:
        print i
        i += 1
        pmid_refs[pmid] = []
        conn.putrequest("GET", "/content/abstract/pubmed_id/"+pmid)
        conn.putheader("Accept", RESPONSE_FORMAT)    
        conn.putheader("X-ELS-APIKey", API_KEY)    
        conn.endheaders()
        response = conn.getresponse()
        root = ET.Element(ET.fromstring(response.read()))
        if response.status != 200:
            paper_errors[pmid] = (response.status, response.reason)
        else:
            for scopus_id in root.tag.iter("refd-itemidlist"):
                status, payload = scopusID_to_pmID(scopus_id[0].text)
                if status:
                    pmid_refs[pmid].append(payload)
                else:
                    ref_errors[scopus_id] = payload                                  
    return pmid_refs

if __name__ == "__main__":
    
    gaf_file = "GOA_Files/gene_association.goa_mouse"
    pmid_list = get_pmids_from_gaf(gaf_file)
    pmid_refs = get_ref_list(pmid_list)
    fHandler = open("SCOPUS_API_Files/"+"pmid_pmid_mouse", 'wb')
    cPickle.dump(pmid_refs, fHandler)
    fHandler.close()
    
    print "Finished Mouse!"
    
    gaf_file = "GOA_Files/gene_association.goa_human"
    pmid_list = get_pmids_from_gaf(gaf_file)
    file = open ("mouse_totla_pmids", 'wb')
    file.write(pmid_list)
    file.close()
#     pmid_refs = get_ref_list(pmid_list)
#     fHandler = open("SCOPUS_API_Files/"+"pmid_pmid_human", 'wb')
#     cPickle.dump(pmid_refs, fHandler)
#     fHandler.close()
    