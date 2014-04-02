#!/usr/bin/env python
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from Bio import Entrez as ez
from Bio.UniProt import GOA
import sys
import os
import argparse
import json
import progressbar

# Logging
from log import Logger

ez.email= "jomaao@miamioh.edu"

pmid_pmid = {}

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
CHROMEDRIVER = os.path.join(CURR_PATH, "Utilities/chromedriver")

CACHED_RESULTS = {}
PAGE_TIMEOUT = 30

def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contained.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    
    pmid_go = {}
    with open(gaf_file, 'r') as unigoa_file:
        pmids = {}
        pmid_prot = {}
        for inrec in GOA.gafiterator(unigoa_file):
            for dbref in inrec['DB:Reference']:
                if dbref[:4] == 'PMID':
                    pmid = dbref[5:]
                    pmids[pmid] = None
                    if pmid not in pmid_go:
                        pmid_go[pmid] = [inrec['GO_ID']]
                    elif inrec['GO_ID'] not in pmid_go[pmid]:
                        pmid_go[pmid].append(inrec['GO_ID'])
                    if pmid not in pmid_prot:
                        pmid_prot[pmid] = [inrec['DB_Object_ID']]
                    elif inrec['DB_Object_ID'] not in pmid_prot[pmid]:
                        pmid_prot[pmid].append(inrec['DB_Object_ID'])
        
    return list(pmids.keys()), pmid_go, pmid_prot



def pmid2doi(pmid_list):
    """
        Divide the pmid list to 200 sublists and query PubMed database for the dois
        @param pmid_list: list of pmids to get the dois for.
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
        @param pmid_list: 200-element list of pmids.
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

def firefox_setup():
    
    fp = webdriver.FirefoxProfile()
    fp.add_extension(extension='Utilities/firebug-1.11.0.xpi')
    fp.set_preference("extensions.firebug.currentVersion", "1.11.0") #Avoid startup screen
    return fp

def get_references(pmids_dois, species, browsername, verbose=0):
    '''
        Starts the process of getting references for a list of dois.        
        @param pmids_dois: Dictionary of pmids and their corresponding dois to search Scopus for their references.
    '''
    
    pmid_refs = {}
    
    if verbose > 0:
        print "Starting selected browser..."

    if browsername.lower() == "chrome":
        os.environ["webdriver.chrome.driver"] = CHROMEDRIVER
        browser = webdriver.Chrome(CHROMEDRIVER)
    elif browsername.lower() == "firefox":
        browser = webdriver.Firefox(firefox_profile=firefox_setup())
    else:
        return Exception('Please choose either FIREFOX or CHROME as your browser')

    count, total = 0, len(pmid_dois.keys())

    bar = progressbar.ProgressBar(maxval=total, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    for pmid, doi in pmids_dois.iteritems():
        logger = Logger(pmid + ".log")
        if pmid not in CACHED_RESULTS:
            pmid_refs[pmid] = parse_scopus_output(search(browser, doi, logger))
            CACHED_RESULTS[pmid] = pmid_refs[pmid]
        else:
            logger.log("Using cached result for this pmid...")
            pmid_refs[pmid] = CACHED_RESULTS[pmid]
        count += 1
        bar.update(count)

    browser.quit()

    bar.finish()

    if verbose > 0:
        print "Writing out results to pmid_pmid_" + species + "_scopus.json..."
    with open("pmid_pmid_" + species + "_scopus.json", 'wb') as out:
        json.dump(pmid_refs, out, indent=4)


def search(browser, doi, logger):
    '''
        Searches one DOI on Scopus and returns a list of references.        
        @param doi: The DOI to search.
    '''
    
    try:
        logger.log("Entering Scopus at url: " + SCOPUS_QUERY_URL)
        browser.get(SCOPUS_QUERY_URL)

        logger.log("Searching for DOI: " + doi)
        search_field = browser.find_element_by_id("searchfield")
        search_field.send_keys(DOI_FORMAT.format(doi))
        
        search_button = browser.find_element_by_class_name("btnSearch")
        search_button.click()
        
        logger.log("Selecting all results from the query... Waiting on page load")
        selectAll = WebDriverWait(browser, PAGE_TIMEOUT).until(EC.presence_of_element_located((By.ID, "allPageCheckBox")))
        selectAll.click()
        
        more_button = browser.find_element_by_id("moreLink")
        more_button.click()
        
        reference_viewer_anchor = browser.find_element_by_css_selector('#moreOptions > div > a.viewReferences')
        reference_viewer_anchor.click()
        

        logger.log("Selecting all the options from the references... Waiting on page load")
        selectAll = browser.find_element_by_id("allPageCheckBox")
        selectAll.click()
        
        export_button = browser.find_element_by_id("export_results")
        export_button.click()

        logger.log("Waiting for export popup...")
         
        format_select = WebDriverWait(browser, PAGE_TIMEOUT).until(EC.presence_of_element_located((By.ID, "custom-radioBtn-Text")))
        format_select.click()

        logger.log("Specifying PMID option only")
        output_select_dropdown = browser.find_element_by_id('exportViewSelectBoxItContainer')
        output_select_dropdown.click()
        
        output_select = browser.find_element_by_css_selector('#exportViewSelectBoxItOptions li:last-child')
        output_select.click()

        uncheck = WebDriverWait(browser, PAGE_TIMEOUT).until(EC.presence_of_element_located((By.ID, "selectedCitationInformationItemsAll-Export")))
        if uncheck.is_selected():
	    uncheck.click()
        
        check = browser.find_element_by_id("selectedBibliographicalInformationItems-Export4")
        if not check.is_selected():
	    check.click()

        logger.log("Triggering the export process")
        export_button = browser.find_element_by_id("exportTrigger")
        export_button.click()
        
        logger.log("Switching browser windows")
        browser.switch_to_window(browser.window_handles[-1])

        pre = WebDriverWait(browser, 200).until(EC.presence_of_element_located((By.TAG_NAME, "pre")))
        text = pre.text.encode('utf-8')
        logger.log(text)
        return text
    except Exception as e:
        logger.log(str(e), "ERROR")
	return ""


def parse_scopus_output(scopus_output):
    '''
        Parses the output of the scopus text format.
        @param scopus_output: the whole text string that Scopus spits out
    '''
    
    references = []
    lines = scopus_output.split("\n")
    for line in lines:
        if line.startswith("PUBMED ID:"):
            references.append(line[11:len(line)-1])  
    return references

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mine scopus for references to papers")   
    parser.add_argument('goa_file', nargs=1, help='The GOA file to use for searching')
    parser.add_argument('-b', '--browser', default='chrome', choices=['chrome', 'firefox'], help='The browser to use for mining')
    parser.add_argument('-s', '--species', help='Name of the species which is being mined')    
    parser.add_argument('-v', '--verbose', action='count', help='Increases the verbosity level. Prints progress information')
    parser.add_argument('-c', '--cache', help="A json file to start from mapping pmids to references")

    args = parser.parse_args()
    goa_file = args.goa_file[0]
    species = args.species if args.species else os.path.splitext(goa_file)[1][1:]

    if args.verbose > 0:
        print "Extracting pmids from the goa file..."
    pmids, pmid_go, pmid_prot = pmids_from_gaf(goa_file)
    
    if args.verbose > 0:
        print "Converting pmids to DOIs for searching..."

    if args.cache and os.path.isfile(args.cache):
        with open(args.cache) as cachedFile:
            # for cached results these are often on a restart so we must
            # remove results with an empty list assuming they might be flawed
            CACHED_RESULTS = { pmid : refs for pmid, refs in json.load(cachedFile).iteritems() if refs }
            
    # don't search for pmids that we will be hitting the cache on
    pmids = [ pmid for pmid in pmids if pmid not in CACHED_RESULTS ]

    pmid_dois = pmid2doi(pmids)
    pmid_dois.update(CACHED_RESULTS)
   
    if args.verbose > 0:
       print "Getting the references for each pmid from Scopus (this could take a long while)..."
 
    get_references(pmid_dois, species, args.browser, args.verbose)
    
