import httplib
import sys

if __name__ == "__main__":
    pmid = sys.argv[1]
    conn = httplib.HTTPConnection("api.elsevier.com:80")
     
    conn.putrequest("GET", "/content/abstract/pubmed_id/"+pmid)
    #conn.putrequest("GET", "/content/search/scopus?query=eid%282-s2.0-0020818830%29+OR+eid%282-s2.0-0033001756%29&field=pubmed-id")
    #conn.putrequest("GET", "/content/search/scopus?query=scopus-id%280030722135%29+OR+scopus-id%2879953118839%29&field=pubmed-id")
     
    #conn.putheader("Accept", "application/json")
    conn.putheader("Accept", "text/xml")
     
    conn.putheader("X-ELS-APIKey", "eeeeac63bdcbd8551deacd2d4d445a00")
     
    conn.endheaders()
    response = conn.getresponse()
    data = response.read()
    file = open("xmlfile.xml", 'wb')
    file.write(data)
    file.close()
    
