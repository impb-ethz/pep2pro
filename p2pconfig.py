#!/usr/bin/env python
"""
p2pconfig.py: pep2pro helper class to read configuration files, 
              validates configuration,
              extract parameters from configuration, 
              checks availability of searchdatabase,
              loads searchdatabase,
              creates locus
"""
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Development"

import configparser
import sys,os
import re
import json

configfile='mascot_search.ini'

def read_config():
	#parse configuration
	config = configparser.ConfigParser()
	if not os.path.isfile(configfile):
		sys.exit("Missing config-file mascot_search.ini.")
	else:
		config.read(configfile)
		return config
		
def extract_regex(config,dbconfig):
	#read and extract regex for loci and decoy 
	#including length
	regextract={}
	regextract['loci']=[e.strip() for e in config.get(dbconfig,'locus').split(',')]
	regextract['locisubstr']=[e.strip() for e in config.get(dbconfig,'locus_len').split(',')]
	regextract['decoy']=[e.strip() for e in config.get(dbconfig,'decoy').split(',')]
	for i in range(0,len(regextract['decoy'])):
		if regextract['decoy'][i]=='True':
			regextract['decoy'][i]=True
		elif regextract['decoy'][i]=='False':
			regextract['decoy'][i]=False
		else:
			sys.exit("Invalid configuration parameter in decoy section:"+regextract['decoy'][i])
	return regextract

def chk_modification(config,mod):
	#check if modification exists in config file
	if not mod in config['MODIFICATIONS']:
		sys.exit("missing modification \""+mod+"\" in config file")
	#initialize dictionary
	ma_mods={}
	#split definition in config 
	tmp_mod=config['MODIFICATIONS'][mod].split(',')
	#first value=amino acid, 2.value=aa with modification,3.value=modification type
	for ii in range(0,len(tmp_mod),3):
		ma_mods[tmp_mod[ii]]=[tmp_mod[ii+1],tmp_mod[ii+2]]
	return ma_mods
		
def chk_searchdb(config,searchdb,chkjson): #chkjson=True during import, False during makep2pdb.py
	#check if fasta search database exists
	dbconfig=''
	#check all databases in configfile
	for db in config['SEARCHDATABASES']:
		if config['SEARCHDATABASES'][db]==searchdb:
			dbconfig=db
			break
	if not dbconfig:
		#the fasta-file is missing in config file
		sys.exit("missing "+searchdb+" in config file")
	if chkjson==True:
		if not os.path.isfile(config['SEARCHDATABASEPATH']['PATH']+"/"+searchdb+".json"):
			#the search
			print "missing index file for "+searchdb+"."
			print "run makep2pdb.py "
			sys.exit(1)
	if chkjson==True:
		if not os.path.isfile(config['SEARCHDATABASEPATH']['PATH']+"/"+searchdb+".stats.json"):
			#the search
			print "missing statistics file for "+searchdb+"."
			print "run makep2pdb.py "
			sys.exit(1)
	return dbconfig	

def load_proteome(config,searchdb):
	#reads and returns created proteome with makep2pdb.py
	proteome={}
	with open(config['SEARCHDATABASEPATH']['PATH']+"/"+searchdb+".json") as data_file:    
		proteome = json.load(data_file)	
	return proteome
	
def makelocus(regexes,protident):
	#inspects protein identifier and returns 
	# 1. locus
	# 2. protein-ident
	# 3. reverse hit flag
	locus={'locus':None,'protident':protident,'decoy':None}
	#is check for valid locus necessary?
	if not regexes['loci'][0]==False:
		#check locus
		for i in range(0,len(regexes['loci'])):
			loc = protident[0:int(regexes['locisubstr'][i])]
			#print i,loc
			regex = re.compile(regexes['loci'][i])
			m=regex.search(loc)
			if m is not None:
				locus['locus']=loc
				locus['decoy']=regexes['decoy'][i]
	return locus
	