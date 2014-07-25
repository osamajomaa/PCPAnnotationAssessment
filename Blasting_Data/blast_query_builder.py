'''
    This file is to build a blast query in a fasta format. The query records are proteins resulting from the 
    intersection of Human genbank files and the Human papers cited by mouse papers. 
'''

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import cPickle

CURR_PATH = "/home/osama/Documents/Research/Project/PCPAnnotationAssessment/"

def make_header(feat_list):
    if not feat_list:
        return "$"
    header = ""
    for feat in feat_list:
        print feat
        header += feat + ";"
    print '----------'
    return header[:-1]

def create_FASTA_file(sprot_file):
    
    handle = open(fname)
    record = SwissProt.parse(handle);
    frecords = []
    for rec in record:
        prot_id = rec.entry_name
        prot_acc = rec.accessions
        refs = rec.references
        prot_refs = []
        for ref in refs:
            for r in ref.references:
                if r[0] == "PubMed":
                    prot_refs.append(r[1])                
        sequence = rec.sequence
        desc = make_header(prot_acc)+"|"+make_header(prot_refs)
        frec = SeqRecord(Seq(sequence, IUPAC.protein), id=prot_id, description=desc)
        frecords.append(frec)
    output_handle = open("human.fasta", "w")
    SeqIO.write(frecords, output_handle, "fasta")
    output_handle.close()
    return "human.fasta"


def uniquify_list(crosspec_rec):
    
    visited = []
    new_crosspec_rec = []
    for pair in crosspec_rec:
        if pair['pmid'] not in visited:
            new_crosspec_rec.append(pair)
            visited.append(pair['pmid'])
    return new_crosspec_rec


def get_intersection(crosspec_pmids, pmid_prot):
    
    inter_prot = {}
    numOfPmids = 0
    numOfProts = 0
    for pmid in crosspec_pmid:
        uniq_crosspec_pmid = uniquify_list(crosspec_pmid[pmid])
        for pair in uniq_crosspec_pmid:
            if pair['species'] == 'Human' or pair['species'] == 'Humans,Mice':
                if pair['pmid'] in pmid_prot:
                    numOfPmids += 1
                    for prot in pmid_prot[pair['pmid']]:
                        #print pair['pmid']
                        if prot not in inter_prot:
                            numOfProts += 1
                            inter_prot[prot] = (pmid_prot[pair['pmid']][prot], [pair['pmid']])
                        elif pair['pmid'] not in inter_prot[prot][1]:
                            inter_prot[prot][1].append(pair['pmid'])
            
    stats_file = open("stats.txt", 'w')
    stats = "Number of papers = " + str(numOfPmids) + "\n"
    stats += "Number of proteins = " + str(numOfProts) + "\n"
    stats_file.write(stats)
    stats_file.close()
    frecords = []
    for prot in inter_prot:
        pmids = inter_prot[prot][1]
        seq = inter_prot[prot][0]
        frec = SeqRecord(Seq(seq, IUPAC.protein), id=prot, description=make_header(pmids))
        frecords.append(frec)
    output_handle = open("query.fasta", "w")
    SeqIO.write(frecords, output_handle, "fasta")
    output_handle.close()
    return "query.fasta"
                           
                           
if __name__ == "__main__":
    crosspec_pmid = cPickle.load(open(os.path.join(CURR_PATH,"mesh_headings/research/mouse_citations")))
    pmid_prot = cPickle.load(open(os.path.join(CURR_PATH,"GO_Proteins/Genbank/HomoSapiens/pmid_prot")))
    get_intersection(crosspec_pmid, pmid_prot);
    
    #query_file = "query.fasta"
    #db_file = "sprot_db/humanblast.db"
    #evalue = 10**-8
    #make_BLAST_query(query_file, db_file, evalue)
    #print "Finish"
    
    
    