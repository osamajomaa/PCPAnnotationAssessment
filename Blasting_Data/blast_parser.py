from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys
import cPickle
import os

CURR_PATH = "/home/jomaao/Project/"
def get_blast_hits(input, hits_output):
    
    bhits = {}
    result_handle = open(os.path.join(CURR_PATH,"Data/blast_results", input))
    records = NCBIXML.parse(result_handle)
    for rec in records:
        for alignment in rec.alignments:
            protein = alignment.title.split('|')[2].split()[1]
            bhits[protein] = None
    hits = list(bhits.keys())

    fHandler = open(os.path.join(CURR_PATH,"Data/BLAST/Hits",hits_output), 'w')
    cPickle.dump(hits, fHandler)
    fHandler.close()
    return hits


def get_blast_misses(bhits, fasta_file, miss_output):
    bmisses = []
    hits = 0
    misses = 0
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        if seq_record.id not in bhits:  
            bmisses.append(seq_record.id)
            misses += 1
        else:
            hits += 1
    
    print "hits = ", hits
    print "misses = ", misses
    
    fHandler = open(os.path.join(CURR_PATH,"Data/BLAST/Misses",miss_output), 'w')
    cPickle.dump(bmisses, fHandler)
    fHandler.close()
    print "Done with misses!"

if __name__ == "__main__":

    input = sys.argv[1]
    ho = sys.argv[2]
    mo = sys.argv[3]
    human_fasta = os.path.join(CURR_PATH,"Data/Swissprot_Files/FASTA_Files/human.fasta")
    #hits = get_blast_hits(input, ho)
    hits = cPickle.load(open(os.path.join(CURR_PATH,"Data/BLAST/Hits",ho)))
    get_blast_misses(hits, human_fasta, mo)
    print "Done!"
    