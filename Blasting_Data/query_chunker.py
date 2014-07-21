'''
    This file is to break down a fasta file into a number of subfiles so that it can be processed 
    in parallel
'''
from Bio import SeqIO
import os

CURR_PATH = "/home/osama/Documents/Research/Project/PCPAnnotationAssessment/Blasting_Data"

def get_size(records):
    size = 0
    for rec in records:
        size += 1
    return size

def fasta_chunker(n, fname):
    records = list(SeqIO.parse(fname, "fasta"))
    length = len(records)
    csize = length/n
    rem = length%n
    start = 0
    end = csize + rem
    i = 1
    while (end <= length):
        chunk = records[start:end]
        output_handle = open(os.path.join(CURR_PATH,"query_chunks/chunk"+str(i)+".fasta"), "w")
        SeqIO.write(chunk, output_handle, "fasta")
        output_handle.close()
        start = end
        end = start + csize
        i += 1

if __name__ == "__main__":
    fastafile = os.path.join(CURR_PATH,"CD-HIT/nr_query.fasta") 
    chunks = 6
    fasta_chunker(chunks, fastafile)
    print "Finish!"
    


    
    