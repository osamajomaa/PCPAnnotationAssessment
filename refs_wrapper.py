import cPickle
import sys
import os

CURR_PATH = os.path.dirname(os.path.realpath(__file__))

def wrap_data(species):
    
    files = os.listdir(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+species+"/Results/data"))
    
    pmid_pmid = {}
    for fname in files:
        par_pmid_pmid = cPickle.load(open(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+species+"/Results/data/"+fname)))
        pmid_pmid.update(par_pmid_pmid)
    
    fHandler = open("SCOPUS_API_Files/"+"pmid_pmid_"+species, 'wb')
    cPickle.dump(pmid_pmid, fHandler)
    fHandler.close()
    
def wrap_errors(species):
    
    files = os.listdir(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+species+"/Results/errors"))
    
    err_data = ""
    err_num = 0
    for fname in files:
        if not fname.endswith("~"):   
            with open(os.path.join(CURR_PATH,"SCOPUS_API_Files/"+species+"/Results/errors/"+fname)) as par_err_file:
                contents = par_err_file.read()
                err_data += contents
                err_num += contents.count("\n")
    
    err_data = "Number of papers not found in Scopus = " + str(err_num) + "\n\n" + err_data
    
    with open(os.path.join(CURR_PATH,"Erronous_Results/not_found_scopus_pmid_"+species+".txt"), 'wb') as err_file:
            err_file.write(err_data)

if __name__ == "__main__":
    
    species = sys.argv[1]
    wrap_data(species)
    wrap_errors(species)
    
    
