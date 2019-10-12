# GET AIRFOIL .DAT FILES FROM UIUC AIRFOIL SITE
# Written by: JoshTheEngineer
# Permission: Dr. Michael Selig (01/16/19)
# Started   : 01/14/19
# Updated   : 01/14/19 - Works as expected
#
# UIUC Airfoil Database: https://m-selig.ae.illinois.edu/ads/coord_database.html

# Importing
from bs4 import BeautifulSoup																	# Import the BeautifulSoup library
import re																						# Import regular expressions
try:																							# Import urllib
    import urllib.request as urllib2
except ImportError:
    import urllib2

# Base filepath for the UIUC airfoil website (used for accessing .dat files)
baseFlpth = "https://m-selig.ae.illinois.edu/ads/coord_seligFmt/"								# Base filepath for saving

# Open the webpage and create the soup
html_page = urllib2.urlopen(baseFlpth)															# Open the URL		
soup      = BeautifulSoup(html_page,'lxml')														# Create the soup

# Loop over all relevant files and save each one
ind   = 1																						# Iteration counter
#links = []																						# Initialize list of links for appending
for link in soup.find_all('a',attrs={'href': re.compile('\.dat', re.IGNORECASE)}):				# Loop over all appropriate links on webpage
    #links.append(link.get('href'))																# Append the link to the list
    
    urllib2.urlretrieve(baseFlpth+link.get('href'), link.get('href').rsplit('/',1)[-1])			# Get the data from the webpage, and save it to the save data file as the link name
    print("Saving file %i" %ind)																# Indicate the link that we are currently saving
    ind = ind + 1																				# Increment the counter
