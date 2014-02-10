#!/usr/bin/env python
from selenium import webdriver
from Bio import Entrez as ez
from Bio.UniProt import GOA
import networkx as nx
import cPickle
import os
import sys

# Graphic stuff
import matplotlib.pyplot as plt

ez.email = "jomaao@miamioh.edu" 


"""
    Global Variables
        pmid_go: Dictionary of all papers in the Uniprot_GOA file and the GO terms in each one.
        pmid_pmid: Dictionary of all papers in the Uniprot_GOA and their references.
        prot_prot: The Protein-Protein Network.
        pcp: The Protein-Citation-Protein Network.
"""

pmid_pmid = {}

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'


CHROMEDRIVER = os.path.join(os.path.dirname(os.path.realpath(__file__)), "chromedriver")

def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contain.
        @ param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    
    pmid_go = {}
    unigoa_file = open(gaf_file)
    pmids = {}
    pmid_prot = {}
    go_terms = []
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                pmids[pmid] = None
                go_terms.append(inrec['GO_ID'])
                if pmid not in pmid_go:
                    pmid_go[pmid] = [inrec['GO_ID']]
                else:
                    pmid_go[pmid].append(inrec['GO_ID'])
                if pmid not in pmid_prot:
                    pmid_prot[pmid] = [inrec['DB_Object_ID']]
                elif inrec['DB_Object_ID'] not in pmid_prot[pmid]:
                    pmid_prot[pmid].append(inrec['DB_Object_ID'])
    # I enforced the list cast here because the dict_key is not subscriptableds))
    return list(pmids.keys()), pmid_go, go_terms, pmid_prot


def remove_high_throughput_papers(pmid_go, pmid_prot):
    
    for pmid in pmid_prot.keys():
        if (len(pmid_prot[pmid]) >= 50):
            del pmid_go[pmid]
            del pmid_prot[pmid]

def pmid2doi(pmid_list):
    """
        Divide the pmid list to 200 sublists and query PubMed database for the dois
        @ param pmid_list: list of pmids to get the dois for.
    """
    
    i = 0 
    j = 200
    pmid_doi = {}
    while (j <= len(pmid_list)):
        pmid_doi = dict(list(pmid_doi.items()) + list(pmid2doi_Helper(pmid_list[i:j]).items()))
        i = j
        j += 200
    j = len(pmid_list)
    pmid_doi = dict(list(pmid_doi.items()) + list(pmid2doi_Helper(pmid_list[i:j]).items()))
    return pmid_doi

def pmid2doi_Helper(pmid_list):
    """
        Query PubMed database to get doi for each pmid which are in 200 sublists
        @ param pmid_list: 200-element list of pmids.
    """
    
    pmid_doi = {}
    pmid = ''
    doi = ''
    pmid_str_list = str(pmid_list)[1:-1]
    handle = ez.efetch("pubmed", id=pmid_str_list, retmode="xml")
    for pubmed_rec in ez.parse(handle):
        for art_id in pubmed_rec['PubmedData']['ArticleIdList']:
            if art_id.attributes['IdType'] == 'doi':
                doi = str(art_id)
            elif art_id.attributes['IdType'] == 'pubmed':
                pmid = str(art_id)        
        if doi != '':
            pmid_doi[pmid] = doi
            doi = ''
    return pmid_doi


class QueryError(Exception):
    pass

def get_references(pmids_dois):
    '''
        Starts the process of getting references for a list of dois.        
        @param pmids_dois: Dictionary of pmids and their corresponding dois to search Scopus for their references.
    '''
    
    pmid_refs = {}
    for pmid, doi in pmids_dois.iteritems():
        pmid_refs[pmid] = parse_scopus_output(search(doi))
    pmid_pmid = pmid_refs
    fHandler = open("pmid_pmid", 'wb')
    cPickle.dump(pmid_pmid, fHandler)
    fHandler.close()
    

def search(doi):
    '''
        Searches one DOI on Scopus and returns a list of references.        
        @param doi: The DOI to search.
    '''
    
    try:
        os.environ["webdriver.chrome.driver"] = CHROMEDRIVER
        browser = webdriver.Chrome(CHROMEDRIVER)
        browser.get(SCOPUS_QUERY_URL)
        search_field = browser.find_element_by_id("searchfield")
        search_field.send_keys(DOI_FORMAT.format(doi))
        
        search_button = browser.find_element_by_class_name("searchButton")
        search_button.click()
        
        selectAll = browser.find_element_by_id("selectAllTop")
        selectAll.click()
        
        more_button = browser.find_element_by_class_name("arrowDownMore")
        more_button.click()
        
        reference_viewer_anchor = browser.find_element_by_xpath('//*[@id="reqMenuList"]/li[1]/a')
        reference_viewer_anchor.click()
        
        selectAll = browser.find_element_by_id("selectAllTop")
        selectAll.click()
        
        export_button = browser.find_element_by_id("export")
        export_button.click()
        
        format_select = browser.find_element_by_css_selector("#exportFormat > option[value=TEXT]")
        format_select.click()
        
        output_select = browser.find_element_by_css_selector("select[name=view] > option[value=SpecifyFields]")
        output_select.click()
        
        uncheck = browser.find_element_by_id("selectedCitationInformationItemsAll-Export")
        uncheck.click()
        
        check = browser.find_element_by_id("selectedBibliographicalInformationItems-Export4")
        check.click()
        
        export_button = browser.find_element_by_css_selector("input[value=Export]")
        export_button.click()
        
        pre = browser.find_element_by_css_selector("pre")
        text = pre.text
        
        browser.quit()
        
        return text
    except:
        browser.quit()
        return ""


def parse_scopus_output(scopus_output):
    '''
        Parses the output of the scopus text format.
    '''
    
    references = []
    outputFile = open("temp_file", 'w')
    outputFile.write(scopus_output)
    outputFile.close()
    outputFile = open("temp_file")
    lines = outputFile.readlines()
    for line in lines:
        if line.startswith("PUBMED ID:"):
            references.append(line[11:len(line)-1])  
    return references


def create_pcp_network(pmid_pmid,pmid_go):
    """
        Create a One Depth Citation Relationship (1DCR) Protein-Citation-Protein network where nodes are GO terms and edges 
        are the hidden relationships between two proteins that are in two papers where one cites the other.
        
    """
    pcp = nx.Graph()
    remove_term = False
    for pmid,go_terms in pmid_go.items():
        for term in go_terms:
            if remove_term == True: # Keep the 1DCR
                go_terms.remove(term)
                remove_term = False
            if pmid_pmid.has_key(pmid):
                for ref_pmid in pmid_pmid[pmid]:
                    if pmid_go.has_key(ref_pmid):
                        for ref_term in pmid_go[ref_pmid]:
                            if not pcp.has_edge(term, ref_term):
                                pcp.add_edge(term, ref_term, co_ocrnce=1)
                                pmid_go[ref_pmid].remove(ref_term) # Keep the 1DCR
                                
                                remove_term = True
                            else:
                                pcp[term][ref_term]['co_ocrnce'] += 1
                                    
    
    return pcp

def create_prot_prot_network(pmid_go):
    """
        Create a Protein-Protein network where nodes are GO terms and edges 
        are the hidden relationships between two proteins that are in the same paper.
    """    
    prot_prot = nx.Graph()
    
    for pmid in pmid_go:
        n = len(pmid_go[pmid])
        for i in range(0,n-1):
            for j in range(i+1,n):
                if not prot_prot.has_edge(pmid_go[pmid][i], pmid_go[pmid][j]):
                    prot_prot.add_edge(pmid_go[pmid][i], pmid_go[pmid][j], co_ocrnce=1)
                else:
                    prot_prot[pmid_go[pmid][i]][pmid_go[pmid][j]]['co_ocrnce'] += 1
    return prot_prot

def get_stats(pmid_pmid, pmid_go, go_terms):
    
    #Number of go terms in dicty file
    num_go_terms = len(go_terms)
    
    #Number of pmids in dicty file
    num_pmids = len(pmid_pmid)
    
    #Avg number of go terms for each pmid
    sum = 0
    for pmid in pmid_go:
        sum += len(pmid_go[pmid])
    avg_go_terms = sum/len(pmid_pmid)
    
    #Avg number of references for each pmid
    sum = 0
    for pmid in pmid_pmid:
        sum += len(pmid_pmid[pmid])
    avg_pmid_ref = sum/len(pmid_pmid)
    
    #Avg number of references for each pmid that are in the dicty file
    sum = 0
    for pmid in pmid_pmid:
        for ref in pmid_pmid[pmid]:
            if ref in pmid_pmid.keys():
                sum += 1
    avg_pmid_dicty_ref = sum/len(pmid_pmid)
    
    #Avg number of occurences in a pmid in the dicty file for each go term
    sum = 0
    for term in go_terms:
        for pmid in pmid_go:
            if term in pmid_go[pmid]:
                sum += 1
    avg_terms_pmid = sum/len(go_terms) 
    
    print "Number of go terms in dicty file = " + str(num_go_terms)
    print "Number of pmids in dicty file = " + str(num_pmids)
    print "Avg number of go terms for each pmid = " + str(avg_go_terms)
    print "Avg number of references for each pmid = " + str(avg_pmid_ref)
    print "Avg number of references for each pmid that are in the dicty file = " + str(avg_pmid_dicty_ref)
    print "Avg number of occurences for each go term in each pmid in dicty file = " + str(avg_terms_pmid)
    
    
    
if __name__ == "__main__":
    
    pmids,pmid_go, go_terms, pmid_prot = pmids_from_gaf("gene_association.goa_dicty")
    #pmid_dois = pmid2doi(pmids)  
    #get_references(pmid_dois)
    
    remove_high_throughput_papers(pmid_go, pmid_prot)
    
    pmid_pmid = cPickle.load(open("/home/jomaao/GitHub/PCPAnnotationAssessment/pmid_pmid"))
    
    prot_prot = create_prot_prot_network(pmid_go)
    pos = nx.graphviz_layout(prot_prot, prog='neato', args='')
    nx.draw(prot_prot,pos,node_size=10,alpha=0.5,node_color="red", with_labels=False)
    plt.savefig("/home/jomaao/GitHub/PCPAnnotationAssessment/prot_prot.png")
    
    pcp = create_pcp_network(pmid_pmid,pmid_go)
    pos = nx.graphviz_layout(pcp, prog='neato', args='')
    nx.draw(pcp,pos,node_size=10,alpha=0.5,node_color="red", with_labels=False)
    plt.savefig("/home/jomaao/GitHub/PCPAnnotationAssessment/prot_citation_prot.png")
    
    get_stats(pmid_pmid, pmid_go, go_terms)
