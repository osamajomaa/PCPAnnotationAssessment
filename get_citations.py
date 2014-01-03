#!/usr/bin/env python
import sys
from Bio import Entrez as ez
from Bio.UniProt import GOA 

ez.email = "jomaao@miamioh.edu" #I changed the email so that if I screwed up, I get the "gentle" message to my email :)
def get_citations(pmid):
    """
    Returns the pmids of the papers this paper cites
    """
    cites_list = []
    handle = ez.efetch("pubmed", id=pmid, retmode="xml")
    pubmed_rec = ez.parse(handle).next()
    if not pubmed_rec['MedlineCitation'].has_key('CommentsCorrectionsList'):
        return []
    for ref in pubmed_rec['MedlineCitation']['CommentsCorrectionsList']:
        if ref.attributes['RefType'] == 'Cites':
            cites_list.append(str(ref['PMID']))
    return cites_list

def pmids_from_gaf(unigoa_file):
    pmids = {}
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                pmids[pmid] = None
    return list(pmids.keys()) # I enforced the list cast here because the dict_key is not subscriptable

def pmid2doi(pmid_list):
    i = 0 
    j = 200
    pmid_doi = {}
    while (j <= len(pmid_list)):
        pmid_doi = dict(list(pmid_doi.items()) + list(pmid2doi_Helper(pmid_list[i:j]).items()))
        i = j
        j += 200
    j = len(pmid_doi)
    pmid_doi = dict(list(pmid_doi.items()) + list(pmid2doi_Helper(pmid_list[i:j]).items()))
    return pmid_doi

def pmid2doi_Helper(pmid_list):
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
        pmid_doi[pmid] = doi
    return pmid_doi
    
def make_scopus_query(doi_list):
    scopus_query = ''
    for doi in doi_list:
        scopus_query += ' DOI("%s") OR \n' % doi
    scopus_query = scopus_query[:-4]
    return scopus_query

def multi_get_citations(pmid_list):
    cites = {}
    i = 0
    for pmid in pmid_list:
        i += 1
        citations = get_citations(pmid)
        if citations == []:
            print i, "poop", pmid
        else:
            cites[pmid] = citations
        if i == 15: break
    return cites
def save_citations(cites,outhandle):
    for citing in cites:
        outhandle.write("%s\t" % citing)
        for pmid in cites[citing]:
            outhandle.write("%s\t" % pmid)
        outfile.write("\n")



if __name__ == '__main__':
    z = get_citations(sys.argv[1])
    for i in z:
        print i
