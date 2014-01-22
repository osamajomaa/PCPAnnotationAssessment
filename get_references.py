#!/usr/bin/env python
from selenium import webdriver
import os

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'

# Chrome driver is os specific download at:
#     https://code.google.com/p/selenium/wiki/ChromeDriver
CHROMEDRIVER = os.path.join(os.path.dirname(os.path.realpath(__file__)), "chromedriver")

class QueryError(Exception):
    pass

def get_references(dois):
    '''
        Starts the process of getting references for a list of dois.
        
        @param dois: List of DOI's to search in Scopus.
    '''
    if not isinstance(dois, list):
        dois = [dois]
    return { doi : parse_scopus_output(search(doi)) for doi in dois }

def search(doi):
    '''
        Searches one DOI on Scopus and returns a list of references.
        
        @param doi: The DOI to search.
    '''
    
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

def parse_scopus_output(scopus_output):
    '''
        Parses the output of the scopus text format.
    '''
    
    return scopus_output

if __name__ == "__main__":
    print get_references("10.1111/j.1432-0436.2003.07109002.x")
    
    
    
    