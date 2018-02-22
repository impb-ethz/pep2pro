# pep2pro
easy pep2pro - from peptide to proteome V2

## Prerequistis
Following python libraries are required:
- json
- click
- configparser
- urllib2
- cookielib
- msparser
- Bio
If you want to import the data into a MySQL Database:
- mysql.connector

## How to use

### setup ini-files
Two ini-files are required:

#### mascot_server.ini
The file contains parameter to access the mascot server (url) and the scripts to call for login and to download files. Additionally it contains the necessary login credentials for the mascot server. 
The file is needed by script get_ma_datfile.py

#### mascot_search.ini
The file contains parameter to filter downloaded mascot dat-files. 

### scripts



