import os
import cPickle

CURR_PATH = "/home/jomaao/Project/"

def wrapper(foldername):
    
    grand_list = []
    files = os.listdir(os.path.join(CURR_PATH, "Data/BLAST", foldername))
    grand_list
    for fname in files:
        sublist = cPickle.load(open(os.path.join(CURR_PATH, "Data/BLAST", foldername, fname)))
        grand_list.extend(sublist)
    
    fhandler = open(os.path.join(CURR_PATH, "Data/BLAST", foldername, "total"+foldername), 'w')
    cPickle.dump(grand_list, fhandler)
    fhandler.close()
    
    print "total number of ", foldername, " = ", len(grand_list)


if __name__ == "__main__":
    wrapper("Hits")
    wrapper("Misses")
    print "Finish!!"
    