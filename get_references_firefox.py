#!/usr/bin/env python
from collections import OrderedDict
from Bio import SwissProt
from selenium import webdriver
from Bio import Entrez as ez
from Bio.UniProt import GOA
import networkx as nx
import operator
import cPickle
import pylab
import os


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
    

def firefox_setup():
    
    fp = webdriver.FirefoxProfile()
    fp.add_extension(extension='Utilities/firebug-1.11.0.xpi')
    fp.set_preference("extensions.firebug.currentVersion", "1.11.0") #Avoid startup screen
    return fp

def search(doi):
    '''
        Searches one DOI on Scopus and returns a list of references.        
        @param doi: The DOI to search.
    '''
    try:
        #os.environ["webdriver.chrome.driver"] = CHROMEDRIVER
        browser = webdriver.Firefox(firefox_profile=firefox_setup())
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

def create_ee_network(pmid_ent):
    """
        Create a Entity-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in the same paper.
        @param pmid_ent: the dictionary of the papers and their entities (in our case protein IDs or GO terms)
    """    
    ee = nx.Graph()
    #ee.add_nodes_from(get_entities(pmid_ent))
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

    
    
    
if __name__ == "__main__":
    
    pmids, pmid_go, pmid_prot = pmids_from_gaf(os.path.join(CURR_PATH,"GOA_Files/gene_association.goa_human"))
    
    pmid_dois = pmid2doi(pmids)  
    
    get_references(pmid_dois)
    
    '''
    print len(pmid_prot)
    remove_high_throughput_papers(pmid_go, pmid_prot,60)
    print len(pmid_prot)
    pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"Pickled_Data/pmid_pmid_dicty")))
    
    
    #GO Terms
    gg_net = create_ee_network(pmid_go)
    draw_network(gg_net, 'gg', os.path.join(CURR_PATH,"gg_net.png"))
    
    gcg_net = create_ece_network(pmid_pmid, pmid_go)
    draw_network(gcg_net, 'gcg', os.path.join(CURR_PATH,"gcg_net.png"))
    
    net_deg_gcg = get_network_degrees(gcg_net)
    draw_network_degrees(OrderedDict(net_deg_gcg.items()[:10]), os.path.join(CURR_PATH,"net_deg_top_gcg.png"), 'Go Terms', "GCG", False)
    draw_network_degrees(net_deg_gcg, os.path.join(CURR_PATH,"net_deg_all_gcg.png"), 'Go Terms', "GCG")
    
    net_deg_gg = get_network_degrees(gg_net)
    draw_network_degrees(OrderedDict(net_deg_gg.items()[:10]), os.path.join(CURR_PATH,"net_deg_top_gg.png"), 'Go Terms', "GG", False)
    draw_network_degrees(net_deg_gg, os.path.join(CURR_PATH,"net_deg_all_gg.png"), 'Go Terms', "GG")
    
    gg_nohigh_net = remove_high_degree_nodes(net_deg_gg, gg_net, 10)
    gcg_nohigh_net = remove_high_degree_nodes(net_deg_gcg, gcg_net, 10)
    draw_network(gg_nohigh_net, 'gg', os.path.join(CURR_PATH,"gg_nohigh_net.png"))
    draw_network(gcg_nohigh_net, 'gcg', os.path.join(CURR_PATH,"gcg_nohigh_net.png"))
    
    
    
    #Poteins
    pp_net = create_ee_network(pmid_prot)
    draw_network(pp_net, 'pp', os.path.join(CURR_PATH,"pp_net.png"))
    
    pcp_net = create_ece_network(pmid_pmid, pmid_prot)
    draw_network(pp_net, 'pcp', os.path.join(CURR_PATH,"pcp_net.png"))
    
    net_deg_pcp = get_network_degrees(pcp_net)
    draw_network_degrees(OrderedDict(net_deg_pcp.items()[:10]), os.path.join(CURR_PATH,"net_deg_top_pcp.png"), 'Proteins', "PCP", False)
    draw_network_degrees(net_deg_pcp, os.path.join(CURR_PATH,"net_deg_all_pcp.png"), 'Proteins', "PCP")
    
    net_deg_pp = get_network_degrees(pp_net)
    draw_network_degrees(OrderedDict(net_deg_pp.items()[:10]), os.path.join(CURR_PATH,"net_deg_top_pp.png"), 'Proteins', "PP", False)
    draw_network_degrees(net_deg_pp, os.path.join(CURR_PATH,"net_deg_all_pp.png"), 'Proteins', "PP")
    
    pp_nohigh_net = remove_high_degree_nodes(net_deg_pp, pp_net, 10)
    pcp_nohigh_net = remove_high_degree_nodes(net_deg_pcp, pcp_net, 10)
    draw_network(pp_nohigh_net, 'pp', os.path.join(CURR_PATH,"pp_nohigh_net.png"))
    draw_network(pcp_nohigh_net, 'pcp', os.path.join(CURR_PATH,"pcp_nohigh_net.png"))
    

'''