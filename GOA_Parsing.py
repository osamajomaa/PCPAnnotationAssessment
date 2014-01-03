import Bio.UniProt.GOA as goa
import sys
import Bio.Entrez as ez

"""""
Retrieving protein references from the yeast association file in GAF 2.0 format
according to different criteria
"""


""""
Retrieve all references cited to annotate proteins with Experimental Evidence Codes 
"""
handle = open("gene_association.goa_yeast")  # open the association gene file of the yeast
proteins = goa.gafiterator(handle) # read all records in the file 
Evidences = {"Evidence":set(["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"])}
print ("GO-annotated proteins supported by Experimental Evidence Code")
for protein in proteins:
    if goa.record_has(protein, Evidences):
        print(protein['DB:Reference'])

""""
Retrieve all references cited to annotate proteins with Experimental Evidence Codes
in the Molecular Function aspect of GO
"""
handle = open("gene_association.goa_yeast")
proteins = goa.gafiterator(handle) 
Evi_Aspect = {"Evidence":set(["EXP", "IDA", "IPI", "IMP", "IGI", "IEP"]), "Aspect":set(["F"])}
print ("GO-annotated proteins supported by Experimental Evidence Code in the Molecular Function Ontology")
for protein in proteins:
    if goa.record_has(protein, Evi_Aspect):
        print(protein['DB:Reference'])
print("\n")

""""
Retrieve all references cited to support the annotation of three UniProt 
proteins whose IDs are: O13516, YPR010C-A, STE2
"""
handle = open("gene_association.goa_yeast")
proteins = list(goa.gafbyproteiniterator(handle)) # get all records in yeast file groupped by object id
IDs = {"DB_Object_ID":set(["A2P2R3","P00045","P50101"])}
for protein in proteins:
    if goa.record_has(protein[0], IDs):
        print("Protein ID:" + protein[0]["DB_Object_ID"]+"\n"+"References:")
        for p in protein:
            print(p['DB:Reference'])
        print("\n")