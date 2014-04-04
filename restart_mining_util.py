import os, sys
import datetime
from time import localtime
import argparse
import json

def filter_log_before_date(files, date, logs_dir='logs'):
    # first find only the log files
    files = [f for f in files if os.path.splitext(f)[1] == '.log']
    d = datetime.datetime.strptime(date, "%B %d %Y %H:%M").timetuple() # Month Day FullYear Hour:Minute
    results = set()
    for f in files:
        filetimesecs = os.path.getmtime(os.path.join(logs_dir, f))
	ftime = localtime(filetimesecs)
	if ftime < d:
	    results.add(os.path.splitext(f)[0])

    return results

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scans the logs for entries before a given date/time' + \
        ' and cross-references with a list of all pmids from a chche file on restart')

    parser.add_argument('cache_file', nargs=1, help='The cache file outputted by closing the browser during the script')
    parser.add_argument('date', nargs=1, help='The date from which to scan the logs in the format: Month Day FullYear')
    parser.add_argument('-l', '--logs', default='logs', help='The log directory to search')
    parser.add_argument('-o', '--outfile', help='Set the outfile for the pmid list. Default is <cache_file_prefix>.restart.json')

    args = parser.parse_args()

    outfile = os.path.splitext(args.cache_file[0])[0] + '.restart.json' if not args.outfile else args.outfile

    with open(args.cache_file[0], 'r') as cache_handle:
        cache = json.load(cache_handle)

    print "Scanning logs directory this could take some time..."

    files = os.listdir(args.logs)
    finished_pmids = filter_log_before_date(files, args.date[0], args.logs)
    pmids = [pmid for pmid, refs in cache.iteritems() if pmid not in finished_pmids]
    
    with open(outfile, 'w') as out:
        json.dump(pmids, out)

