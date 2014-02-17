#!/usr/bin/env python
from selenium import webdriver
from Bio import Entrez as ez
from Bio.UniProt import GOA
import networkx as nx
import cPickle
import os
import operator
import pylab
from collections import OrderedDict

# Graphic stuff
import matplotlib.pyplot as plt

ez.email = "jomaao@miamioh.edu" 

pmid_pmid = {}

SCOPUS_QUERY_URL = "http://www.scopus.com/search/form.url?display=advanced&clear=t&origin=searchbasic"
DOI_FORMAT = 'DOI("{0}")'

CURR_PATH = os.path.dirname(os.path.realpath(__file__))
CHROMEDRIVER = os.path.join(CURR_PATH, "chromedriver")


def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contain.
        @ param unigoa_file: Uniprot_GOA association file in gaf format. 
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
    # I enforced the list cast here because the dict_key is not subscriptableds))
        
    return list(pmids.keys()), pmid_go, pmid_prot


def remove_high_throughput_papers(pmid_go, pmid_prot, threshold):

    for pmid in pmid_prot.keys():
        if (len(pmid_prot[pmid]) >= threshold):
            del pmid_go[pmid]
            del pmid_prot[pmid]
    

def get_network_degrees(network):
       
    ent_deg = OrderedDict()
    
    for node in network.nodes():
        ent_deg[node] = len(network.neighbors(node))
    
    ent_deg = OrderedDict(sorted(ent_deg.iteritems(), key=operator.itemgetter(1), reverse=True))
    
    return ent_deg

def draw_network_degrees(net_deg, path_name, xtitle, net_name, all_nodes=True):
    
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




def create_ece_network(pmid_pmid,pmid_ent):
    """
        Create a One Depth Citation Relationship (1DCR) Entity-Citation-Entity network where nodes are entities and edges 
        are the hidden relationships between two entities that are in two papers where one cites the other.
        
    """
    
    ece = nx.Graph()
    ece.add_nodes_from(get_entities(pmid_ent))
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
    """    
    ee = nx.Graph()
    ee.add_nodes_from(get_entities(pmid_ent))
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
    
    top_deg = OrderedDict(hd_nodes.items()[:topx])
    for node in top_deg:
        network.remove_node(node)
    
    return network

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
    
    pmids,pmid_go, pmid_prot = pmids_from_gaf(os.path.join(CURR_PATH,"gene_association.goa_dicty"))
    #pmid_dois = pmid2doi(pmids)  
    #get_references(pmid_dois)
    
    print len(pmid_prot)
    remove_high_throughput_papers(pmid_go, pmid_prot,100)
    print len(pmid_prot)
    pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"pmid_pmid")))
    
    
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
    

