# pep2pro
easy pep2pro - from peptide to proteome V2

## Prerequistis
The scripts were written for **Python Version 2**.

Following python libraries are required:
* json
* click
* configparser
* urllib2
* cookielib
* msparser (s. http://www.matrixscience.com/msparser_download.html)
* Bio

If you want to import the data into a MySQL Database:

* mysql.connector

## Installation

    $ git clone https://github.com/impb-ethz/pep2pro.git


## How to use

### setup ini-files
Two ini-files are required:

#### mascot_server.ini
The file contains parameter to access the mascot server (url) and the scripts to call for login and to download files. Additionally it contains the necessary login credentials for the mascot server. 
The file is needed by script get_ma_datfile.py
```
[SERVER]
url=http://<url-to-mascot-server>/mascot/
login=cgi/login.pl
file=x-cgi/ms-status.exe

[USER]
username=<mascot username>
password=<mascot password>

[DATFILEPATH]
;path to which datfiles will be downloaded  
path=/var/mascot/datfiles/
```

#### mascot_search.ini
The file contains parameter to parse and load mascot dat-files. 
```
[PEPTIDESUMMARY]
;defines the minimal probability of the ms_peptidesummary object call
MinProb=0.1
;defines the maxHitsToReport proteins, whatever the protein score. 0 means no cutoff
MaxHitsToReport=0 

[SCANTITLE]
;create regular expression to extract measurement filename from query.title
;e.g.: 426: Scan 2248 (rt=41.4859) [S:\datafolder\20100421_05_sample4.RAW]
v1 = \[.\:(.*)\]
;e.g.: 219: Scan 2899 (rt=24.279) [\\file-server.ethz.ch\Data2San\20140829_elute_no.raw]                	
v2 = \[\\\\file\-server\.ethzh\.ch\\Data2San(.*)\]

[MODIFICATIONS]
;define modifications as used in dat-file
;with amino acid and new value and modification type
;separated by comma, e.g. Phospho (S)=S,S166,Phospho.
;if n-term modification, use . (dot) so modification 
;will be added at the beginning of modified peptide
;IMPORTANT!!!!!!!!!!!!
;if modification type is specified, the algorithm
;checks if peptides (unmodified) between rank 1 and rank 2 hits
;are the same. If that is the case and the ratio of ionsscore is <= 0.4
;no modification is shown as peptide_modified because of
;uncertain position, but the peptide is marked as modified 
;in the specified types, e.g. Phospho
;IMPORTANT!!!!!!!!!!!!
;for N-term modification, only one triplet is allowed, 
;e.g. Acetyl (Protein N-term)=.,42,Acetyl
;for others, multiple triplets can be defined, 
;e.g. Phospho_STY (STY)=S,S166,Phospho,T,T181,Phospho,Y,Y243,Phospho
Acetyl (Protein N-term)=.,42,Acetyl
Phospho_STY (STY)=S,S166,Phospho,T,T181,Phospho,Y,Y243,Phospho

[SEARCHDATABASES]
;define all search database with an alias.
;the alias has to exists as SECTION in the ini file
;to define regular expression to determine valid 
;locus and if protein is decoy
;shortversion should not be longer then 25 characters
tair10=TAIR10_20110117.fasta 
                                                                           
[SEARCHDATABASEPATH]                                                       
;path where the search databases used in mascot can be found              
;and where the summary files created with make2pdb.py will be stored      
path=/var/databases/searchdbs                                             
                                                                          
[MYSQLDB]                                                                  
;connection parameters and credentials to store data in                    
;a MySQL Database                                                          
user=<database username>                                                  
password=<database password>                                               
host=<mysql server>                                                        
database=<database name>                                                  
                                                                          
[tair10]
;define all valid loci as regular expression, comma separated
;1. regex for the locus, 2. substring of the locus, 3. decoy loci False/True
;e.g. AT1G01010.1 is a valid locus
;REV_AT1G01010.1 is a valid locus and decoy
;zz_RAV03.881.12 is not a valid locus
locus=^AT.G.....$,^REV_AT.G.....$
locus_len=9,13
decoy=False,True
```

### How To

1. **makep2pdb.py**

Run makep2pdb.py to generate needed summary files for the database you used in mascot. This creates two additional files *inputfile*.json and *inputfile*.stats.json which will be saved in the [SEARCHDATABASEPATH] (s. mascot_search.ini).

  $ python make2pdb.py <searchdatabase>

2. **get_ma_datfile.py**

This script downloads a .dat file from the mascot server for further processing with parse_ma_datfile.py. The server url and access credentials have to be provided in mascot_server.ini.
```sh
python get_ma_datfile.py <mascot-directory> <mascot-datfile>
```
3. **parse_ma_datfile.py**

This script filters a .dat file according to minimum ionsscore, maximum peptide expectation value and the pep2pro ambiguity filter. The script creates an outputfile in json format as *local-datfile*.json
```sh
python parse_ma_datfile.py <local-datfile> <min-ionsscore> <max-expectation-value>
```
4. **load_ma_datfile.py**

This script imports the json resultfile from parse_ma_datfile.py into the MySQL Database, specified in the MYSQLDB section of the mascot_search.ini.
```sh
python load_ma_datfile.py <local-json-resultfile>
```
Optional parameters are 
* experiment
* genotype
* tissue
* treatment
* replicate
* fraction

if not specified during script initiation, e.g. --experiement mytestexperiment --genotype col0 etc., you will be prompted for the input.

## License
The msparser library is subject to the licese agreement of Mascot, s. Mascot Parser Licence Agreement at http://www.matrixscience.com/msparser_download.html 

Copyright &copy; 2018, Matthias Hirsch-Hoffmann, ETH Zurich. Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html) or see [LICENSE.txt](LICENSE.txt).


## Disclamer
This software is supplied 'as is' without any warranty or guarantee of support. ETH Zurich is not responsible for its use, misuse, or functionality. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability arising from, out of, or in connection with this software.
