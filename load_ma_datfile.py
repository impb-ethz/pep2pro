#!/usr/bin/env python
"""
load_ma_datfile.py: pep2pro script to load analyzed mascot datfile into mysql database.
"""
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Development"

import sys
import time
import os.path
import p2pconfig
from p2pmysql import pep2prodb
import re
import json
#argumentparser
import click

config=p2pconfig.read_config()

@click.command()
#which arguments and options are needed? 
#arguments=mandatory, option=possible
#*************** ARGUMENETS *****************************
#input-file
@click.argument('ma_datfile', type=click.Path(exists=True,readable=True))
@click.option('--experiment', prompt='Please enter an experiment name',help='Experiment name.',default='')
@click.option('--genotype'  , prompt='Please enter a sample genotype ',help='Sample genotype.',default='')
@click.option('--tissue'    , prompt='Please enter a sample tissue ',help='Sample tissue.',default='')
@click.option('--treatment' , prompt='Please enter a sample treatment ',help='Sample treatment.',default='')
@click.option('--replicate' , prompt='Please enter a replicate ',help='Replicate.',default='')
@click.option('--fraction'  , prompt='Please enter a fraction ',help='Fraction.',default='')
def main(ma_datfile,experiment,genotype,tissue,treatment,replicate,fraction):
	#check json data
	with open(ma_datfile) as data_file:    
		#load data
		data = json.load(data_file)	
		#check search database
		searchdb = p2pconfig.chk_searchdb(config,data['fastaver'],True)
		#get database search id
		searchdbid=p2p.get_searchdatabase_id(searchdb)	
		if searchdbid == None:
			#create search database
			searchdbid=p2p.create_searchdatabase_id(searchdb,data['fastaver'])
			#import search-databaes stats json file
			if p2p.load_searchdatabase(config,searchdbid,data['fastaver']):
				p2p.mysql_commit()
		#get proteins
		p2p.get_proteins(searchdbid)	
		#get datfile_id
		datfile_id=p2p.get_datfile_id(data['datfile'],data['ionsscore'],data['probability'],searchdbid,experiment,genotype,tissue,treatment,replicate,fraction)
		for hit in data['hits']:
			#load hits
			peptide_id,peptide_modified_id=p2p.get_pepetide_peptide_modified_id(hit['goodaccessions'],hit['peptide'],hit['misscleaved'],hit['peptide_modified'],hit['modifications'])
			measurement_id=p2p.update_measurement(hit['measurement'])
			#insert hit record
			sql = ("insert into hit (datfile_id,measurement_id,peptide_id,peptide_modified_id "
					",decoy,ambiguous,multiloci,query,rank,ionsscore,probability,scannbr"
					",charge,retention_time,calcmass,deltamass,obsmass,expmass) "
					" values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s )")
			record=[datfile_id,measurement_id,peptide_id,peptide_modified_id
			,hit['decoyhit'],hit['ambiguous'],hit['multiloci'],hit['query'],hit['rank']
			,hit['ionsscore'],hit['probability'],hit['scannbr']
			,hit['charge'],hit['rettime'],hit['calcmass'],hit['deltamass']
			,hit['obsmass'],hit['expmass']
			]
			csr=p2p.mysql_execute(sql,record,1)
			hit_id=csr.lastrowid
			#add hit modification type
			p2p.update_hit_modificationtypes(hit_id,hit['modificationtypes'])
				
	
if __name__ == '__main__':
	p2p = pep2prodb(config['MYSQLDB']['user']
		,config['MYSQLDB']['password']
		,config['MYSQLDB']['host']
		,config['MYSQLDB']['database'])
	main()
	p2p.close()


