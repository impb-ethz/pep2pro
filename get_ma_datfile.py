#!/usr/bin/env python
"""
get_ma_datfile.py: pep2pro script to download mascot datfile from mascot server.
"""
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Development"

import sys
import click
import configparser
import urllib2
from cookielib import CookieJar

@click.command()
#which arguments and options are needed? 
#arguments=mandatory, option=possible
#*************** ARGUMENETS *****************************
#input-file
@click.argument('ma_dir')
@click.argument('ma_datfile')

def main(ma_dir,ma_datfile):
	#parse configuration
	config = configparser.ConfigParser()
	config.read('mascot_server.ini')
	
	cj = CookieJar()
	opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
	#set necessary cookies
	params="username="+config['USER']['username']+"&action=login&password="+config['USER']['password']
	
	response = opener.open(config['SERVER']['url']+config['SERVER']['login'], params)	

	msaction=config['SERVER']['url']+config['SERVER']['file']+"?Autorefresh=false&Show=RESULTFILE&DateDir="+ma_dir+"&ResJob="+ma_datfile+"&BrowserSafe=true"

	content=opener.open(msaction).read()
	f = open(config['DATFILEPATH']['path']+"/"+ma_datfile, 'w')
	f.write(content)
	f.close()
    	
if __name__ == '__main__':
	main()




