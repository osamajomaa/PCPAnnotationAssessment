#!/usr/bin/env python
import sys
from Bio import Entrez as ez
import Bio.UniProt.GOA as goa


ez.email = "jomaao@miamioh.edu"
def get_citations(pmid):
    """
    Returns the pmids of the papers this paper cites
    """
    cites_list = []
    handle = ez.efetch("pubmed", id=pmid, retmode="xml")
    pubmed_rec = ez.parse(handle).__next__()
    for ref in pubmed_rec['MedlineCitation']['CommentsCorrectionsList']:
        if ref.attributes['RefType'] == 'Cites':
            cites_list.append(str(ref['PMID']))
    return cites_list


f = open ("papers and citations.txt","w")
st = "GO-annotated proteins supported by IGI evidence (Inferred from Genetic Interaction)\n"
handle = open("gene_association.goa_yeast")
proteins = goa.gafiterator(handle) 
Evi_Aspect = {"Evidence":set(["IGI"])}
for protein in proteins:
    if goa.record_has(protein, Evi_Aspect):
        for p in protein['DB:Reference']:
            if p[:4] == "PMID":
                st += "Main PubMed reference: "+ p +"\n"
                citations = get_citations(p[5:])
                for cit in citations:
                    st += cit + "  "
                st += "\n"
f.write(st)
f.close()
        
        
        