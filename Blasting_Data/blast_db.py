from Bio.Blast.Applications import NcbiblastpCommandline
from subprocess import call
import sys
    
E_VALUE = 10**-8

def build_BLAST_db(path, fasta_file, blast_db):
    cmd= "makeblastdb -in " + fasta_file + " -dbtype 'prot' -out " + blast_db
    call(cmd, shell=True)
    
def make_BLAST_query(query_file, db_file, out_file):
    cline = NcbiblastpCommandline(query=query_file, db=db_file, outfmt=5, evalue=E_VALUE, out=out_file)
    print (cline)
    cline()

if __name__ == "__main__":
    
    query_file = sys.argv[1]
    db_file = sys.argv[2]
    out_file = sys.argv[3]
    #query_file = "query.fasta"
    #db_file = "humanblast.db"
    make_BLAST_query(query_file, db_file, out_file)
    print "Finish!!"