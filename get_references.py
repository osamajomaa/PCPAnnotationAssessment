#!/usr/bin/env python
from selenium import webdriver
from Bio import Entrez as ez
from Bio.UniProt import GOA
import networkx as nx
import pickle
import os

ez.email = "jomaao@miamioh.edu" 


"""
    Global Variables
        pmid_go: Dictionary of all papers in the Uniprot_GOA file and the GO terms in each one.
        pmid_pmid: Dictionary of all papers in the Uniprot_GOA and their references.
        PP: The Protein-Protein Network.
        PCP: The Protein-Citation-Protein Network.
"""

pmid_go = {}
pmid_pmid = {}
PP = nx.Graph()
PCP = nx.Graph()

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'


CHROMEDRIVER = os.path.join(os.path.dirname(os.path.realpath(__file__)), "chromedriver")

def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contain.
        @ param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    
    unigoa_file = open(gaf_file)
    pmids = {}
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                pmids[pmid] = None
                if pmid not in pmid_go:
                    pmid_go[pmid] = [inrec['GO_ID'][3:]]
                else:
                    pmid_go[pmid].append(inrec['GO_ID'][3:])
    return list(pmids.keys()) # I enforced the list cast here because the dict_key is not subscriptable

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
    pickle.dump(pmid_pmid, fHandler)
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


def create_PCP_network():
    """
        Create a One Depth Citation Relationship (1DCR) Protein-Citation-Protein network where nodes are GO terms and edges 
        are the hidden relationships between two proteins that are in two papers where one cites the other.
        
    """
    
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
                            if not PCP.has_edge(term, ref_term):
                                PCP.add_edge(term, ref_term, co_ocrnce=1)
                                pmid_go[ref_pmid].remove(ref_term) # Keep the 1DCR
                                remove_term = True
                            else:
                                PCP[term][ref_term]['co_ocrnce'] += 1
                                    
    

def create_PP_network():
    """
        Create a Protein-Protein network where nodes are GO terms and edges 
        are the hidden relationships between two proteins that are in the same paper.
    """    
    
    for pmid in pmid_go:
        n = len(pmid_go[pmid])
        for i in range(0,n-1):
            for j in range(i+1,n):
                if not PP.has_edge(pmid_go[pmid][i], pmid_go[pmid][j]):
                    PP.add_edge(pmid_go[pmid][i], pmid_go[pmid][j], co_ocrnce=1)
                else:
                    PP[pmid_go[pmid][i]][pmid_go[pmid][j]]['co_ocrnce'] += 1
                

if __name__ == "__main__":
    
    pmids = pmids_from_gaf("gene_association.goa_dicty")
    pmid_dois = pmid2doi(pmids)  
    #get_references(pmid_dois)
    d_file = open("pmid_pmid")
    pmid_pmid = pickle.load(d_file)
    create_PP_network()
    create_PCP_network()
    
    