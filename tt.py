import cPickle

pmid_pmid = cPickle.load(open("SCOPUS_API_Files/"+"pmid_pmid_mouse"))

for pmid in pmid_pmid:
    print pmid
    print pmid_pmid[pmid]
    print "\n\n"