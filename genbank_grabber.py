import sys
import cPickle
import io
import os
import unicodedata
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "jomaao@miamioh.edu"
CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def grab_files(organism):
    
    handle = Entrez.esearch(db="protein", term=organism+"[Organism]", retmax=10000000, retmode="xml")
    record = Entrez.read(handle)
    i = 0
    for item_id in record['IdList']:
        try:
            hndl = Entrez.efetch(db="protein", id=item_id, rettype="gb", retmode="text")
        except Exception:
            pass
        else:
            try:
                protein_desc = hndl.read()
            except Exception:
                pass
            else:
                file = open("Genbank_Files/"+organism+"/protein"+str(i)+".gbk", "w")
                file.write(protein_desc)
                file.close()
                i += 1

def parse_files(organism):
    
    pmid_prot = {}
    protein_files = os.listdir(os.path.join(CURR_PATH,"Genbank_Files/"+organism))
    for fname in protein_files:        
        fqname = os.path.join(CURR_PATH,"Genbank_Files/"+organism+"/"+fname)
        try:
            protFile = SeqIO.parse(fqname, "genbank")
        except Exception:
            pass
        else:
            try:
                for record in protFile:
                    protein = record.name
                    seq = record.seq._data
                    if protein.strip() != "" and 'references' in record.annotations:
                        for ref in record.annotations['references']:
                            try:
                                pmid = ref.pubmed_id
                            except Exception:
                                pass
                            else:
                                if pmid.strip() != "":
                                    if pmid not in pmid_prot:
                                            pmid_prot[pmid] = {protein:seq}
                                    else:
                                        pmid_prot[pmid][protein] = seq
            except Exception:
                pass
    
    fHandler = open("GO_Proteins/Genbank/"+organism+"/"+"pmid_prot", 'wb')
    cPickle.dump(pmid_prot, fHandler)
    fHandler.close()

if __name__ == "__main__":
    organism = sys.argv[1]
    #grab_files(organism)
    parse_files(organism)
    print "Finish Successfully"
    