#!/usr/bin/env python
from collections import OrderedDict
from go_clustering import load_data
from go_clustering import pickle_data
from selenium import webdriver
from Bio import Entrez as ez
from Bio import SwissProt
from Bio.UniProt import GOA
import networkx as nx
from Bio import Medline
import operator
import cPickle
import pylab
import json
import sys
import os



# Logging
from log import Logger

# Graphic stuff
import matplotlib.pyplot as plt


ez.email = "jomaao@miamioh.edu" 

pmid_pmid = {}

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
CHROMEDRIVER = os.path.join(CURR_PATH, "Utilities/chromedriver")


def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contained.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    
    pmid_go = {}
    unigoa_file = open(gaf_file)
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


def remove_high_throughput_papers(pmid_go, pmid_prot, threshold):
    '''
        Remove the high throughput papers according to a threshold
        @param pmid_go: A dictionary of papers and the GO terms they are annotated to
        @param pmid_prot: A dictionary of papers and the proteins they contain  
    '''    
    for pmid in pmid_prot.keys():
        if (len(pmid_prot[pmid]) >= threshold):
            del pmid_go[pmid]
            del pmid_prot[pmid]
    

def get_network_degrees(network):
    '''
        Get the degrees of nodes in a EE or ECE network
        @param network: The network to get calculate the degrees for its nodes
    '''
    ent_deg = OrderedDict()
    
    for node in network.nodes():
        ent_deg[node] = len(network.neighbors(node))
    
    ent_deg = OrderedDict(sorted(ent_deg.iteritems(), key=operator.itemgetter(1), reverse=True))
    
    return ent_deg

def create_desc_file(net_deg, net_type):
    pass
    '''
    handle = open("uniprot_sprot_19_feb_2014.dat")
    
    net_deg_desc = OrderedDict()
    for prot in net_deg.keys():
        net_deg_desc[prot] = [net_deg[prot],'']
    
    
    desc_text = "{0[0]:<15}{0[1]:<15}{0[2]:<5}".format("Protein Degree Description".split()) + "\n"
    line = [0,0,0] 
    for record in SwissProt.parse(handle):
        for acc in record.accessions:
            if acc in net_deg_desc.keys():
                net_deg_desc[acc][1] = record.description
                break
        
    for prot in net_deg_desc.keys():
        line[0] = prot
        line[1] = net_deg_desc[prot][0]
        line[2] = net_deg_desc[prot][1]
        desc_text += "{0[0]:<15}{0[1]:<15}{0[2]:<5}".format(line) + "\n"
        

    if net_type == "PP":
        desc_file = open("PP.desc", 'w')
    elif net_type == "PCP":
        desc_file = open("PCP.desc", 'w')
        
    desc_file.write(desc_text)
    '''
    
def draw_network_degrees(net_deg, path_name, xtitle, net_name, all_nodes=True):
    '''
        Plot a histogram of the top N nodes having the highest degrees in a network
        @param net_deg: dictionary of nodes in the network with their corresponding degree
        @param path_name: the path to store the exported graph image to
        @param xtitle: the title of the x axis in the graph
        @param net_name: the network name: EE or ECE
        @param all_nodes: a flag to determine whether to draw the whole nodes or just the top N 
    '''
    
    xaxis = list(range(0,len(net_deg)))
    yaxis = list(net_deg.values())
    xvals = list(net_deg.keys())
        
    
    if (all_nodes):
        pylab.plot(xaxis, yaxis)
        pylab.fill_between(xaxis,yaxis,0,color='cyan')
        pylab.xticks(visible = False)
        pylab.title("Degrees of nodes in the " + net_name + " " + "network")
        
        
    else:
        pylab.bar(xaxis, yaxis, align='center')
        pylab.gcf().subplots_adjust(bottom=0.25)
        pylab.xticks(xaxis, xvals, rotation=90)
        pylab.title("Top nodes that have the highest degree in the " + net_name + " " + "network")        
        if net_name in ["PP", "PCP"]:
            create_desc_file(net_deg, net_name)
              
    pylab.xlabel(xtitle)
    pylab.ylabel('Network Degree')        
    pylab.grid(True)
    pylab.tight_layout()
    pylab.savefig(path_name)
    pylab.close()
 
        
        

def get_entities (pmid_ent):
    
    entities = {}
    for pmid in pmid_ent.keys():
        for ent in pmid_ent[pmid]:
            entities[ent] = None
    
    return list(entities.keys())
        


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

def get_references(pmids_dois, species, webBrowser='chrome'):
    '''
        Starts the process of getting references for a list of dois.        
        @param pmids_dois: Dictionary of pmids and their corresponding dois to search Scopus for their references.
    '''
    
    pmid_refs = {}
    
    if webBrowser.lower() == "chrome":
        os.environ["webdriver.chrome.driver"] = CHROMEDRIVER
        browser = webdriver.Chrome(CHROMEDRIVER)
    elif webBrowser.lower() == "firefox":
        browser = webdriver.Firefox(firefox_profile=firefox_setup())
    else:
        return Exception('Please choose either FIREFOX or CHROME as your browser')

    for pmid, doi in pmids_dois.iteritems():
        logger = Logger(pmid + ".log")
        pmid_refs[pmid] = parse_scopus_output(search(browser, doi, logger))

    browser.quit()
    pmid_pmid = pmid_refs
    fHandler = open("pmid_pmid_"+species, 'wb')
    cPickle.dump(pmid_pmid, fHandler)
    fHandler.close()
    

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
        
        search_button = browser.find_element_by_class_name("searchButton")
        search_button.click()
        
        logger.log("Selecting all results from the query...")
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

        logger.log("Moving to export screen")
        
        format_select = browser.find_element_by_css_selector("#exportFormat > option[value=TEXT]")
        format_select.click()
        

        logger.log("Specifying PMID option only")
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
        
        return text
    except Exception as e:
        return ""


def parse_scopus_output(scopus_output):
    '''
        Parses the output of the scopus text format.
        @param scopus_output: the whole text string that Scopus spits out
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


def create_ece_network(pmid_pmid,pmid_ent):
    """
        Create a One Depth Citation Relationship (1DCR) Entity-Citation-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in two papers where one cites the other.
        @param pmid_pmid: the dictionary of the paper and their ciations
        @param pmid_ent: the dictionary of the papers and their entites (in our case protein IDs or GO terms)
    """
    
    ece = nx.Graph()
    #ece.add_nodes_from(get_entities(pmid_ent))
    visited = []
    for pmid in pmid_ent.keys():
        if pmid not in visited:
            for ent in pmid_ent[pmid]:
                if pmid_pmid.has_key(pmid):
                    for ref_pmid in pmid_pmid[pmid]:
                        if pmid_ent.has_key(ref_pmid) and ref_pmid not in visited:        
                            for ref_ent in pmid_ent[ref_pmid]:
                                if ece.has_edge(ent, ref_ent):
                                    ece[ent][ref_ent]['co_ocrnce'] += 1
                                else:
                                    ece.add_edge(ent, ref_ent, co_ocrnce=1)
                            visited.append(ref_pmid)
            visited.append(pmid)

    return ece

def create_mch_network(mpmid_hpmid, mpmid_ent, hpmid_ent):
    """
        Create a One Depth Citation Relationship (1DCR) Entity-Citation-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in two papers where one cites the other.
        @param pmid_pmid: the dictionary of the paper and their ciations
        @param pmid_ent: the dictionary of the papers and their entites (in our case protein IDs or GO terms)
    """
    
    ece = nx.Graph()
    #ece.add_nodes_from(get_entities(pmid_ent))
    visited = []
    for mpmid in mpmid_ent.keys():
        if mpmid not in visited:
            for ment in mpmid_ent[mpmid]:
                if mpmid_hpmid.has_key(mpmid):
                    for hpmid in mpmid_hpmid[mpmid]:
                        if hpmid_ent.has_key(hpmid) and hpmid not in visited:        
                            for hent in hpmid_ent[hpmid]:
                                if ece.has_edge(ment, hent):
                                    ece[ment][hent]['co_ocrnce'] += 1
                                else:
                                    ece.add_edge(ment, hent, co_ocrnce=1)
                            visited.append(hpmid)
            visited.append(mpmid)

    return ece


def create_hcm_network(hpmid_mpmid, hpmid_ent, mpmid_ent):
    """
        Create a One Depth Citation Relationship (1DCR) Entity-Citation-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in two papers where one cites the other.
        @param pmid_pmid: the dictionary of the paper and their ciations
        @param pmid_ent: the dictionary of the papers and their entites (in our case protein IDs or GO terms)
    """
    
    ece = nx.Graph()
    #ece.add_nodes_from(get_entities(pmid_ent))
    visited = []
    for hpmid in hpmid_ent.keys():
        if hpmid not in visited:
            for hent in hpmid_ent[hpmid]:
                if hpmid_mpmid.has_key(hpmid):
                    for mpmid in hpmid_mpmid[hpmid]:
                        if mpmid_ent.has_key(mpmid) and mpmid not in visited:        
                            for ment in mpmid_ent[mpmid]:
                                if ece.has_edge(hent, ment):
                                    ece[hent][ment]['co_ocrnce'] += 1
                                else:
                                    ece.add_edge(hent, ment, co_ocrnce=1)
                            visited.append(mpmid)
            visited.append(hpmid)

    return ece


def create_ee_network(pmid_ent):
    """
        Create a Entity-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in the same paper.
        @param pmid_ent: the dictionary of the papers and their entities (in our case protein IDs or GO terms)
    """    
    ee = nx.Graph()
    #ee.add_nodes_from(get_entitihpmid_entes(pmid_ent))
    for pmid in pmid_ent:
        n = len(pmid_ent[pmid])
        for i in range(0,n-1):
            for j in range(i+1,n):
                if not ee.has_edge(pmid_ent[pmid][i], pmid_ent[pmid][j]):
                    ee.add_edge(pmid_ent[pmid][i], pmid_ent[pmid][j], co_ocrnce=1)
                else:
                    ee[pmid_ent[pmid][i]][pmid_ent[pmid][j]]['co_ocrnce'] += 1
    return ee


def draw_network(network, net_name, image_path):
    '''
        Draw a network
        @param network: the network to draw.
        @param net_name: the name of the network; EE or ECE.
        @param imag_path: the path to store the exported graph image to.
    '''
    
    if net_name == 'pp':
        plt.title("Protein-Protein Network")
    elif net_name == 'pcp':
        plt.title("Protein-Citation-Protein Network")
    elif net_name == 'gg':
        plt.title("GO-GO Network")
    elif net_name == 'gcg':
        plt.title("GO-Citation-GO Network")
    elif net_name == 'gcg_hcm':
        plt.title("GO-Citation-GO Network - Human Citing Mouse")
    elif net_name == 'gcg_mch':
        plt.title("GO-Citation-GO Network - Mouse Citing Human")
    elif net_name == 'pcp_hcm':
        plt.title("Protein-Citation-Protein Network - Human Citing Mouse")
    elif net_name == 'pcp_mch':
        plt.title("Protein-Citation-Protein Network - Mouse Citing Human")
    
    pos = nx.graphviz_layout(network, prog='neato', args='')
    nx.draw(network,pos,node_size=10,alpha=0.5,node_color="red", with_labels=False)
    
    plt.tight_layout()
    plt.savefig(image_path)
    plt.close()
    

def remove_high_degree_nodes(hd_nodes, network, topx):
    '''
        Remove the N highest degree nodes from a network.
        @param hd_nodes: the dictionary of nodes and their corresponding degrees.
        @param network: the network to remove the highest N degree nodes from.
        @param topx: N; the number of nodes to remove from the network.
    '''
    top_deg = OrderedDict(hd_nodes.items()[:topx])
    for node in top_deg:
        network.remove_node(node)
    
    for node in network.nodes():
        if network.degree(node) < 1:
            network.remove_node(node)
    
    return network


def human_mouse_cross_ref(hu_pmid_pmid, mo_pmid_pmid):
    
    hpmid_mpmid = {}
    mpmid_hpmid = {}
    
    for hpmid in hu_pmid_pmid.keys():
        hpmid_mpmid[hpmid] = []
        for ref_hpmid in hu_pmid_pmid[hpmid]:
            if ref_hpmid in mo_pmid_pmid:
                hpmid_mpmid[hpmid].append(ref_hpmid)
    
    for pmid in hpmid_mpmid.keys():
        print pmid + "            " + str(hpmid_mpmid[pmid])
    for mpmid in mo_pmid_pmid.keys():
        mpmid_hpmid[mpmid] = []
        for ref_mpmid in mo_pmid_pmid[mpmid]:
            if ref_mpmid in hu_pmid_pmid:
                mpmid_hpmid[mpmid].append(ref_mpmid)
    
    return hpmid_mpmid, mpmid_hpmid

def write_cross_ref(spec1, spec2, s1pmid_s2pmid, file_path):
    
    ref_num = 0    
    desc_text = "{0[0]:<15}{0[1]:<15}".format([spec1, spec2])
    desc_text += "\n"
    for pmid in s1pmid_s2pmid.keys():
        citations = ""
        for ref_pmid in s1pmid_s2pmid[pmid]:
            citations += ref_pmid + ', '
            ref_num += 1
        desc_text += "{0[0]:<15}{0[1]:<15}".format([pmid, citations[:-2]]) + "\n"
    desc_text = "Number of " + spec1 + " to " + spec2 + " citations = " + str(ref_num) + "\n" + desc_text
    
    desc_file = open(file_path, 'w')
    desc_file.write(desc_text)
            

def to_list(data):
    for key in data.keys():
        data[key] = list(data[key]) 


def search_species(meshHeadings):
    desc_names = []
    for heading in meshHeadings:
        for qual_name in heading['QualifierName']:
            desc_names.append(qual_name)
        desc_names.append(heading['DescriptorName'])
    return desc_names

def search_headings(mesh_hds, term):
    for heading in mesh_hds:
        if heading.lower().find(term) != -1:
            return True
    return False

def add_to_crosspec(crosspec_pmid, pubmed_val, desc_names, pmid):
    if 'Humans' in desc_names and 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans,Mice'})
    elif 'Humans' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Humans'})
    elif 'Mice' in desc_names:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'Mice'})
    else:
        crosspec_pmid[pmid].append({'pmid':pubmed_val, 'species':'None'})
    return crosspec_pmid
        

def medline(pmid_pmid, spec1, spec2):
    
    #crosspec_pmid = {}
    crosspec_pmid_methods = {}
    crosspec_pmid_review = {}
    crosspec_pmid_research = {}
    
    for pmid in pmid_pmid.keys():
        try:
            handle = ez.efetch(db="pubmed", id=pmid_pmid[pmid], rettype="medline", retmode="xml")
            crosspec_pmid_methods[pmid] = []
            crosspec_pmid_review[pmid] = []
            crosspec_pmid_research[pmid] = []
            records = ez.parse(handle)
            
            for pubmed_rec in records:
                desc_names = search_species(pubmed_rec['MedlineCitation']['MeshHeadingList'])
                if search_headings(desc_names, 'methods'):
                    crosspec_pmid_methods = add_to_crosspec(crosspec_pmid_methods, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)
                elif search_headings(desc_names, 'review'):
                    crosspec_pmid_review = add_to_crosspec(crosspec_pmid_review, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)
                else:
                    crosspec_pmid_research = add_to_crosspec(crosspec_pmid_research, pubmed_rec['MedlineCitation']['PMID'], desc_names, pmid)                    
        except Exception:
            pass
    
    pickle_data(crosspec_pmid_methods, os.path.join(CURR_PATH, "mesh_headings/methods/"+spec1+"_citations"))
    pickle_data(crosspec_pmid_review, os.path.join(CURR_PATH, "mesh_headings/review/"+spec1+"_citations"))
    pickle_data(crosspec_pmid_research, os.path.join(CURR_PATH, "mesh_headings/research/"+spec1+"_citations"))


def uniquify_list(crosspec_rec):
    
    visited = []
    new_crosspec_rec = []
    for pair in crosspec_rec:
        if pair['pmid'] not in visited:
            new_crosspec_rec.append(pair)
            visited.append(pair['pmid'])
    return new_crosspec_rec

def build_crosspec_files(crosspec_pmid, spec1, spec2, paper_type):
    desc_text_s1s1 = "{0[0]:<15}{0[1]:<15}".format([spec1, spec1])+ "\n"
    desc_text_s1s2 = "{0[0]:<15}{0[1]:<15}".format([spec1, spec2])+ "\n"
    desc_text_s1n = "{0[0]:<15}{0[1]:<15}".format([spec1, 'Other'])+ "\n"
    count_s1s1 = 0
    count_s1s2 = 0
    count_s1n = 0
    for pmid in crosspec_pmid.keys():
        cit_s1s1 = ""
        cit_s1s2 = ""
        cit_s1n = ""
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['species'] == spec1:
                cit_s1s1 += pair['pmid'] + ", "
                count_s1s1 += 1
            elif pair['species'] == spec2:
                cit_s1s2 += pair['pmid'] + ", "
                count_s1s2 += 1
            elif pair['species'] == 'None':
                cit_s1n += pair['pmid'] + ", "
                count_s1n += 1
            elif pair['species'] == spec1+","+spec2:
                cit_s1s2 += pair['pmid'] + ", " 
                cit_s1s1 += pair['pmid'] + ", "
                count_s1s1 += 1
                count_s1s2 += 1        
        desc_text_s1s1 += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1s1[:-2]]) + "\n"
        desc_text_s1s2 += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1s2[:-2]]) + "\n"
        desc_text_s1n += "{0[0]:<15}{0[1]:<15}".format([pmid, cit_s1n[:-2]]) + "\n"
    desc_text_s1s1 = "Number of " + spec1 + " to " + spec1 + " citations = " + str(count_s1s1) + "\n\n" + desc_text_s1s1
    desc_text_s1s2 = "Number of " + spec1 + " to " + spec2 + " citations = " + str(count_s1s2) + "\n\n" + desc_text_s1s2
    desc_text_s1n  = "Number of " + spec1 + " to " + "Other" + " citations = " + str(count_s1n) + "\n\n" + desc_text_s1n
    
    desc_file_s1s1 = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_"+spec1), 'w')
    desc_file_s1s1.write(desc_text_s1s1)
    desc_file_s1s2 = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_"+spec2), 'w')
    desc_file_s1s2.write(desc_text_s1s2)
    desc_file_s1n = open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/"+spec1+"/"+spec1+"_cit_Other"), 'w')
    desc_file_s1n.write(desc_text_s1n)

def draw_histograms(crosspec_pmid, spec1, spec2, paper_type):
    
    s1s1_data, s1s2_data, s1n_data = get_histogram_data(crosspec_pmid, spec1, spec2)
    
    draw_histogram_helper(s1s1_data, spec1, spec1, 'red', paper_type)
    draw_histogram_helper(s1s2_data, spec1, spec2, 'green', paper_type)
    draw_histogram_helper(s1n_data, spec1, 'Others', 'blue', paper_type)

def draw_cited_histograms(crosspec_pmid, spec1, spec2, top, paper_type):
    
    s1s1_data, s1s2_data, s1n_data = get_cited_histogram_data(crosspec_pmid, spec1, spec2)
    
    draw_cited_histogram_helper(s1s1_data, spec1, spec1, top, 'red', paper_type)
    draw_cited_histogram_helper(s1s2_data, spec1, spec2, top, 'green', paper_type)
    draw_cited_histogram_helper(s1n_data, spec1, 'Others', top, 'blue', paper_type)
    
    
def draw_histogram_helper(hist_data, spec1, spec2, bcolor, paper_type):
    
    xaxis = []
    yaxis = []
    for key,value in hist_data.items():
        yaxis.append(value[0])
        xaxis.append(key)
    
    #yaxis = hist_data.values()[0]
    
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
    

def get_histogram_data(crosspec_pmid, spec1, spec2):
    
    s1s1_hist_data = {}
    s1s2_hist_data = {}
    s1n_hist_data = {}
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
                
    s1s1_hist_data = dict(sorted(s1s1_hist_data.items()))
    s1s2_hist_data = dict(sorted(s1s2_hist_data.items()))
    s1n_hist_data = dict(sorted(s1n_hist_data.items()))
    return s1s1_hist_data, s1s2_hist_data, s1n_hist_data 

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
    
        
    
if __name__ == "__main__":
    '''
    hu_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_pmid_human.json")))
    mo_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/mouse/pmid_pmid_mouse.json")))

    medline(mo_pmid_pmid, 'mouse', 'human')
    medline(hu_pmid_pmid, 'human', 'mouse')
    '''
    
    paper_type = 'research'
    spec2 = 'Mice'
    spec1 = 'Humans'
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"mesh_headings/"+paper_type+"/human_citations")))
    build_crosspec_files(crosspec_pmid, spec1, spec2, paper_type)
    draw_cited_histograms(crosspec_pmid, spec1, spec2, 10, paper_type)
    draw_histograms(crosspec_pmid, spec1, spec2, paper_type)
    
    
    '''
    #Mouse Citing Human      
    
    onto = 'all'
    species = 'M_Cit_H'
    hu_pmid_go, hu_pmid_prot, hu_go_prot = load_data('human', onto=onto)
    mo_pmid_go, mo_pmid_prot, mo_go_prot = load_data('mouse', onto=onto)    
    remove_high_throughput_papers(hu_pmid_go, hu_pmid_prot,60)
    remove_high_throughput_papers(mo_pmid_go, mo_pmid_prot, 60)
    hu_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_pmid_human.json")))
    mo_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/mouse/pmid_pmid_mouse.json")))
    hpmid_mpmid, mpmid_hpmid = human_mouse_cross_ref(hu_pmid_pmid, mo_pmid_pmid)
    '''

    
    '''
    #GO Terms
    
    go_mch_net = create_mch_network(mpmid_hpmid, mo_pmid_go, hu_pmid_go)
    draw_network(go_mch_net, 'gcg_mch', os.path.join(CURR_PATH,species+"/go_mch_net.png"))
    
    net_deg_go_mch = get_network_degrees(go_mch_net)
    draw_network_degrees(OrderedDict(net_deg_go_mch.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_go_mch.png"), 'Go Terms', "GCG - Mouse Citing Human", False)
    draw_network_degrees(net_deg_go_mch, os.path.join(CURR_PATH,species+"/net_deg_all_go_mch.png"), 'Go Terms', "GCG - Mouse Citing Human")
    
    go_mch_nohigh_net = remove_high_degree_nodes(net_deg_go_mch, go_mch_net, 10)
    draw_network(go_mch_nohigh_net, 'gcg_mch', os.path.join(CURR_PATH,species+"/go_mch_nohigh_net.png"))
    
    
    
    #Poteins
    p_mch_net = create_mch_network(mpmid_hpmid, mo_pmid_prot, hu_pmid_prot)
    draw_network(p_mch_net, 'pcp_mch', os.path.join(CURR_PATH,species+"/p_mch_net.png"))
    
    net_deg_p_mch = get_network_degrees(p_mch_net)
    draw_network_degrees(OrderedDict(net_deg_p_mch.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_p_mch.png"), 'Proteins', "PCP - Mouse Citing Human", False)
    draw_network_degrees(net_deg_p_mch, os.path.join(CURR_PATH,species+"/net_deg_all_p_mch.png"), 'Proteins', "PCP - Mouse Citing Human")
    
    
    p_mch_nohigh_net = remove_high_degree_nodes(net_deg_p_mch, p_mch_net, 10)
    draw_network(p_mch_nohigh_net, 'pcp_mch', os.path.join(CURR_PATH,species+"/p_mch_nohigh_net.png"))
    
    '''
    '''
    #Human Citing Mouse
    onto = 'all'
    species = 'H_Cit_M'
    hu_pmid_go, hu_pmid_prot, hu_go_prot = load_data('human', onto=onto)
    mo_pmid_go, mo_pmid_prot, mo_go_prot = load_data('mouse', onto=onto)
    remove_high_throughput_papers(hu_pmid_go, hu_pmid_prot,60)
    remove_high_throughput_papers(mo_pmid_go, mo_pmid_prot, 60)
    hu_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_pmid_human.json")))
    mo_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/mouse/pmid_pmid_mouse.json")))
    hpmid_mpmid, mpmid_hpmid = human_mouse_cross_ref(hu_pmid_pmid, mo_pmid_pmid)
    
    #GO Terms
    
    go_hcm_net = create_hcm_network(hpmid_mpmid, hu_pmid_go, mo_pmid_go)
    draw_network(go_hcm_net, 'gcg_hcm', os.path.join(CURR_PATH,species+"/go_hcm_net.png"))
    
    net_deg_go_hcm = get_network_degrees(go_hcm_net)
    draw_network_degrees(OrderedDict(net_deg_go_hcm.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_go_hcm.png"), 'Go Terms', "GCG - Human Citing Mouse", False)
    draw_network_degrees(net_deg_go_hcm, os.path.join(CURR_PATH,species+"/net_deg_all_go_hcm.png"), 'Go Terms', "GCG - Human Citing Mouse")
    
    go_hcm_nohigh_net = remove_high_degree_nodes(net_deg_go_hcm, go_hcm_net, 10)
    draw_network(go_hcm_nohigh_net, 'gcg_hcm', os.path.join(CURR_PATH,species+"/go_hcm_nohigh_net.png"))

    
    
    #Poteins
    p_hcm_net = create_hcm_network(hpmid_mpmid, hu_pmid_prot, mo_pmid_prot)
    draw_network(p_hcm_net, 'pcp_hcm', os.path.join(CURR_PATH,species+"/p_hcm_net.png"))
    
    net_deg_p_hcm = get_network_degrees(p_hcm_net)
    draw_network_degrees(OrderedDict(net_deg_p_hcm.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_p_hcm.png"), 'Proteins', "PCP - Human Citing Mouse", False)
    draw_network_degrees(net_deg_p_hcm, os.path.join(CURR_PATH,species+"/net_deg_all_p_hcm.png"), 'Proteins', "PCP - Human Citing Mouse")
    
    
    p_hcm_nohigh_net = remove_high_degree_nodes(net_deg_p_hcm, p_hcm_net, 10)
    draw_network(p_hcm_nohigh_net, 'pcp_hcm', os.path.join(CURR_PATH,species+"/p_hcm_nohigh_net.png"))
    
'''
    
    
    '''
    species = sys.argv[1][21:]
    
    pmids, pmid_go, pmid_prot = pmids_from_gaf(os.path.join(CURR_PATH,"GOA_Files/"+sys.argv[1]))
    
    pmid_dois = pmid2doi(pmids)  
    
    if (len(sys.argv) == 3):
        get_references(pmid_dois, species, sys.argv[2])
    elif (len(sys.argv) == 2):
        get_references(pmid_dois, species)
    
    '''
    
    '''
    onto = 'all'
    species = 'mouse'
    pmid_go, pmid_prot, go_prot = load_data(species, onto=onto)

    to_list(pmid_go)
    to_list(pmid_prot)
    to_list(go_prot)
    
    print len(pmid_prot)
    remove_high_throughput_papers(pmid_go, pmid_prot,60)
    print len(pmid_prot)
    #pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_pmid_"+species+".json")))
    pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/"+species+"/pmid_pmid_"+species+".json")))
    
    #GO Terms
    gg_net = create_ee_network(pmid_go)
    draw_network(gg_net, 'gg', os.path.join(CURR_PATH,species+"/gg_net.png"))
    
    gcg_net = create_ece_network(pmid_pmid, pmid_go)
    draw_network(gcg_net, 'gcg', os.path.join(CURR_PATH,species+"/gcg_net.png"))
    
    net_deg_gcg = get_network_degrees(gcg_net)
    draw_network_degrees(OrderedDict(net_deg_gcg.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_gcg.png"), 'Go Terms', "GCG", False)
    draw_network_degrees(net_deg_gcg, os.path.join(CURR_PATH,species+"/net_deg_all_gcg.png"), 'Go Terms', "GCG")
     
    net_deg_gg = get_network_degrees(gg_net)
    draw_network_degrees(OrderedDict(net_deg_gg.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_gg.png"), 'Go Terms', "GG", False)
    draw_network_degrees(net_deg_gg, os.path.join(CURR_PATH,species+"/net_deg_all_gg.png"), 'Go Terms', "GG")
     
    gg_nohigh_net = remove_high_degree_nodes(net_deg_gg, gg_net, 10)
    gcg_nohigh_net = remove_high_degree_nodes(net_deg_gcg, gcg_net, 10)
    draw_network(gg_nohigh_net, 'gg', os.path.join(CURR_PATH,species+"/gg_nohigh_net.png"))
    draw_network(gcg_nohigh_net, 'gcg', os.path.join(CURR_PATH,species+"/gcg_nohigh_net.png"))

    #Poteins
    pp_net = create_ee_network(pmid_prot)
    draw_network(pp_net, 'pp', os.path.join(CURR_PATH,species+"/pp_net.png"))
    
    pcp_net = create_ece_network(pmid_pmid, pmid_prot)
    draw_network(pcp_net, 'pcp', os.path.join(CURR_PATH,species+"/pcp_net.png"))
    
    net_deg_pcp = get_network_degrees(pcp_net)
    draw_network_degrees(OrderedDict(net_deg_pcp.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_pcp.png"), 'Proteins', "PCP", False)
    draw_network_degrees(net_deg_pcp, os.path.join(CURR_PATH,species+"/net_deg_all_pcp.png"), 'Proteins', "PCP")
    
    net_deg_pp = get_network_degrees(pp_net)
    draw_network_degrees(OrderedDict(net_deg_pp.items()[:10]), os.path.join(CURR_PATH,species+"/net_deg_top_pp.png"), 'Proteins', "PP", False)
    draw_network_degrees(net_deg_pp, os.path.join(CURR_PATH,species+"/net_deg_all_pp.png"), 'Proteins', "PP")
    
    pp_nohigh_net = remove_high_degree_nodes(net_deg_pp, pp_net, 10)
    pcp_nohigh_net = remove_high_degree_nodes(net_deg_pcp, pcp_net, 10)
    draw_network(pp_nohigh_net, 'pp', os.path.join(CURR_PATH,species+"/pp_nohigh_net.png"))
    draw_network(pcp_nohigh_net, 'pcp', os.path.join(CURR_PATH,species+"/pcp_nohigh_net.png"))
    
'''


#     hu_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_pmid_human.json")))
#     mo_pmid_pmid = json.load(open(os.path.join(CURR_PATH,"Pickled_Data/mouse/pmid_pmid_mouse.json")))
#     pmid_go = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/human/pmid_go_human")))
#     
#     total_mouse_pmid = []
#     for pmid in mo_pmid_pmid.keys():
#         total_mouse_pmid.append(pmid)
#         total_mouse_pmid.extend(mo_pmid_pmid[pmid])
#     
#     print len(total_mouse_pmid)
#     print len(mo_pmid_pmid)
#     x = 0
#     for pmid in hu_pmid_pmid.keys():
#         if pmid in mo_pmid_pmid:
#             print pmid
#     print x
#     
#     hpmid_mpmid, mpmid_hpmid = human_mouse_cross_ref(hu_pmid_pmid, mo_pmid_pmid)
#     write_cross_ref('Human', 'Mouse', hpmid_mpmid, os.path.join(CURR_PATH,"H_Cit_M/hu_cit_mo.txt"))
#     write_cross_ref('Mouse', 'Human', mpmid_hpmid, os.path.join(CURR_PATH,"M_Cit_H/mo_cit_hu.txt"))






