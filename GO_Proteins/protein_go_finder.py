from collections import OrderedDict
from Bio.UniProt import GOA
import cPickle
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def mostciting_sort(pmid_cits):
    cit_count = OrderedDict()
    for pmid in pmid_cits:
        cit_num = len(pmid_cits[pmid])
        if cit_num not in cit_count:
            cit_count[cit_num] = [pmid]
        else:
            cit_count[cit_num].append(pmid)
    
    cit_count = OrderedDict(sorted(cit_count.items(), reverse=True))
    for k in cit_count:
        print k
    return cit_count

def pmids_from_gaf(gaf_file):
    """
        Get the papers cited in the Uniprot_GOA file by their PMID and get the GO terms each paper contained.
        @param unigoa_file: Uniprot_GOA association file in gaf format. 
    """
    unigoa_file = open(gaf_file)
    pmid_go = {}
    pmid_prot = {}
    for inrec in GOA.gafiterator(unigoa_file):
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                if pmid not in pmid_go:
                    pmid_go[pmid] = [inrec['GO_ID']]
                elif inrec['GO_ID'] not in pmid_go[pmid]:
                    pmid_go[pmid].append(inrec['GO_ID'])
                if pmid not in pmid_prot:
                    pmid_prot[pmid] = [inrec['DB_Object_ID']]
                elif inrec['DB_Object_ID'] not in pmid_prot[pmid]:
                    pmid_prot[pmid].append(inrec['DB_Object_ID'])
        
    return pmid_go, pmid_prot

def build_goprot_file(pmid_count, pmid_go, pmid_prot, species):
    
    file_contents = ""
    for count in pmid_count:
        for pmid in pmid_count[count]:
            file_contents += "PMID: " + pmid + "\n"
            file_contents += "Number of Citations = " + str(count) + "\n"
            file_contents += "Proteins: \n"
            for prot in pmid_prot[pmid]:
                file_contents += prot + ", "
            if len(pmid_prot[pmid]) > 0:
                file_contents = file_contents[:-2]
            file_contents += "\nGO Terms:\n"
            for go in pmid_go[pmid]:
                file_contents += go + ", "
            if len(pmid_go[pmid]) > 0:
                file_contents = file_contents[:-2]
            file_contents += "\n\n\n\n"
    
    file = open(os.path.join(CURR_PATH,"GO_Proteins/GAF/"+species+"_goprot.txt"), 'w')
    file.write(file_contents)
    file.close()

if __name__ == "__main__":
    spec = 'human'
    pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+"pmid_pmid_"+spec)))
    pmid_count = mostciting_sort(pmid_pmid)
    gaf_file = os.path.join(CURR_PATH,"GOA_Files/gene_association.goa_"+spec)
    pmid_go, pmid_prot = pmids_from_gaf(gaf_file)
    build_goprot_file(pmid_count, pmid_go, pmid_prot, spec)
    
    
    