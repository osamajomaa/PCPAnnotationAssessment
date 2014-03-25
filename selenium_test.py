from selenium import webdriver
import os
import sys

SCOPUS_QUERY_URL = "http://www.google.com"
CURR_PATH = os.path.dirname(os.path.realpath(__file__))
CHROMEDRIVER = os.path.join(CURR_PATH, "Utilities/chromedriver")

def firefox_setup():
    
    fp = webdriver.FirefoxProfile()
    fp.add_extension(extension='Utilities/firebug-1.11.0.xpi')
    fp.set_preference("extensions.firebug.currentVersion", "1.11.0") #Avoid startup screen
    return fp

def search_site(toyText, webBrowser='chrome'):

    if webBrowser.lower() == "chrome":
        os.environ["webdriver.chrome.driver"] = CHROMEDRIVER
        browser = webdriver.Chrome(CHROMEDRIVER)
    elif webBrowser.lower() == "firefox":
        browser = webdriver.Firefox(firefox_profile=firefox_setup())
    else:
        return Exception('Please choose either FIREFOX or CHROME as your browser')


    return search(browser, toyText)

    #browser.quit()

def search(browser, searchTxt):
    
    try:

        browser.get(SCOPUS_QUERY_URL)

        search_field = browser.find_element_by_id("gbqfq")
        search_field.send_keys(searchTxt)
        
        search_button = browser.find_element_by_id("gbqfb")
        search_button.click()
        
        search_result = browser.find_element_by_id("gbqfq")
        text = search_result.get_attribute('value')
        
        return text
    
    except Exception as e:
        return ""

if __name__ == "__main__":
    
    if (len(sys.argv) == 2):
        print search_site(sys.argv[1])
    elif (len(sys.argv) == 3):
        print search_site(sys.argv[1], sys.argv[2])
    else:
        print "Error!!!"
    