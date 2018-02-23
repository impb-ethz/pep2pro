#!/usr/bin/env python
"""
parse_ma_datfile.py: pep2pro script to parse mascot datfile.
"""
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Development"

import msparser
import sys
import time
import os.path
import p2pconfig
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
@click.argument('min_ionsscore',type=int)
@click.argument('max_probability',type=float)
def main(ma_datfile,min_ionsscore,max_probability):
	scanformat=False
	result={}
	result['datfile']=os.path.basename(ma_datfile)
	result['ionsscore']=min_ionsscore
	result['probability']=max_probability
	result['hits']=[]
	file = msparser.ms_mascotresfile(ma_datfile)
	if file.isValid() : 
		if file.isMSMS() :
			#get fasta file and check if data
			result['fastaver']=file.getFastaVer()
			dbconfig=p2pconfig.chk_searchdb(config,result['fastaver'],True)
			
			#load regular expressions for database
			regexes=p2pconfig.extract_regex(config,dbconfig)
			#load proteome for ambiguity filter
			proteome = p2pconfig.load_proteome(config,result['fastaver'])
			
			ma_mods={} # variable to hold modification abbreviation
			params = file.params()
			# read modification to define modification abbreviation from config file
			i = 1
			while params.getVarModsName(i):
				mod=params.getVarModsName(i) 
				ma_mods[i]=p2pconfig.chk_modification(config,mod)
				i+=1 #increase counter		
			# all modifications were found in configuration file
			#print ma_mods
			#sys.exit(1)
			#start reading mascot datfile			
			results=msparser.ms_peptidesummary(file, msparser.ms_mascotresults.MSRES_MAXHITS_OVERRIDES_MINPROB, float(config['PEPTIDESUMMARY']['MinProb']), int(config['PEPTIDESUMMARY']['MaxHitsToReport']))	
#			print file.getNumQueries()
#			print file.getNumHits(msparser.ms_mascotresfile.SEC_SUMMARY)
#			print results.getNumberOfHits()
#			print results.getMaxRankValue()
			for query in range(1,file.getNumQueries()): # all queries
				p2p_query={}
				#print "Query: "+str(query)
				for rank in range(1,results.getMaxRankValue()): #all ranks
					#print "Rank: "+str(rank)
					#get Peptide for query/rank
					pept = results.getPeptide(query,rank)
					#get div. values for peptide
					peptide = pept.getPeptideStr()
					ionsscore = pept.getIonsScore()
					peptide_mod = pept.getVarModsStr()
					#check ionscore
					if peptide and ionsscore > min_ionsscore: 
						#get probability
						probability = results.getPeptideExpectationValue(pept.getIonsScore(),query)
						#check probability
						if probability < max_probability: 
							#check rank 1 
							if rank==1:	
								p2p_query['query']=query
								p2p_query['rank']=rank
								p2p_query['peptide']=peptide
								p2p_query['ionsscore']=ionsscore
								p2p_query['probability']=probability
								
								p2p_query['peptide_mod'] = peptide_mod
								
								p2p_query['charge'] = pept.getCharge()
								
								p2p_query['calcmass']=pept.getMrCalc()
								p2p_query['deltamass']=pept.getDelta()
								p2p_query['obsmass']=pept.getObserved()
								p2p_query['expmass']=pept.getMrExperimental()
								
								iq = msparser.ms_inputquery(file,query)
								
								title=iq.getStringTitle(1)
								
								if not scanformat:
									scanformat=scantitleformat(title)
									#print scanformat
									
								#print query,rank,title,peptide,peptide_mod,ionsscore,probability					
								p2p_query['measurement']=scantitle(title,scanformat)
	
								p2p_query['scannbr'] = iq.getScanNumbers()
								p2p_query['rettime'] = iq.getRetentionTimes()
								numprot=pept.getNumProteins()
								#print scannbr,rettime,numprot
								p2p_query['accessions']=[]
								for i in range(1,(numprot+1)):
									prot = pept.getProtein(i)
									accession = prot.getAccession()
									p2p_query['accessions'].append(p2pconfig.makelocus(regexes,accession))
			# AMBIGUITY FILTER *******************************
								p2p_query['decoyhit']=0 #accession is from decoy 
								p2p_query['ambiguous']=0 #peptide maps to multiple proteins/loci
								p2p_query['multiloci']=0 #peptide maps to multiple loci with same protein
								multiloci=0
								p2p_query['goodaccessions']=[] #list to hold all good accessions 
			# decoy HIT ************************************
			# only hits with single decoy accession
								if len(p2p_query['accessions'])==1:
									if p2p_query['accessions'][0]['locus']==None:
										#next hit as locus is none
										#print "locus None"
										break
									if p2p_query['accessions'][0]['decoy']==True:
										p2p_query['decoyhit']=1
									#add accession to goodaccessions list
									p2p_query['goodaccessions'].append(p2p_query['accessions'][0])
								else:
			# AMBIGUOUS OR MULTILOCI*****************************
									#loop all accessions
									for acc in p2p_query['accessions']:
										#check if decoy result
										if acc['decoy']==False:
											#check if valid locus
											if not acc['locus'] == None:
												#append valid accession to goodaccessions list
												p2p_query['goodaccessions'].append(acc) 
									#check if goodaccessions is has an element
									if len(p2p_query['goodaccessions']) == 0:
										#no goodaccessions locus, or all decoy
										break
									if len(p2p_query['goodaccessions']) > 1:
										#verify all accession with all others
										for i in range(0,len(p2p_query['goodaccessions'])-1):
											for ii in range(i+1,len(p2p_query['goodaccessions'])):
												#check if accessions share protein sequence
												if p2p_query['goodaccessions'][ii]['protident'] in proteome[p2p_query['goodaccessions'][i]['protident']]:
													#count mulitlocis
													multiloci+=1
										#calculate gausssche sum. if identical with multiloci
										#all accessions share the same protein sequence				
										gsum= (len(p2p_query['goodaccessions'])-1) * len(p2p_query['goodaccessions']) / 2 
										if gsum == multiloci:
											#mark peptide as multiloci
											p2p_query['multiloci']=1
										else:
											#mark peptide as ambiguous
											p2p_query['ambiguous']=1
											if multiloci > 0:
												p2p_query['multiloci']=1
								#last dictionary element for check with rank 2 peptide
								p2p_query['removemods']=False
								#print "pass"
							#check rank 2
							elif rank==2: 
								#check rank 1 vs rank 2 hit if peptide sequence are identical
								#and remove modifications if ratio between ionsscore is <= 0.4.
								#print "rank 2 check "+"="*30
								if peptide==p2p_query['peptide']:
									ratio = (ionsscore/p2p_query['ionsscore'])
									if (ionsscore/p2p_query['ionsscore'])<=0.4:
										#remove phosphorylations
										p2p_query['removemods']=True
										#print "remove mods"
							else: #avoid records rank > 2
								#print "break - "*10
								break
						else:
							#print "bad probability:"+str(probability)
							break
					else:
						#print "bad ionsscore:"+str(ionsscore)
						break		
					if rank == 2:
						#break query loop if rank == 2
						#print "rank 2 break"
						break
				if p2p_query:
					if len(p2p_query['goodaccessions']) > 0:
						del p2p_query['accessions'] # not needed anymore
						#create modified peptide, modificationtypes and details on modification
						p2p_query['peptide_modified'],p2p_query['modificationtypes'],p2p_query['modifications'],p2p_query['misscleaved']=make_modpeptide(ma_mods,p2p_query['peptide'],p2p_query['peptide_mod'],p2p_query['removemods'])
						del p2p_query['peptide_mod'] # not needed anymore
						del p2p_query['removemods'] # not needed anymore
						result['hits'].append(p2p_query)
			#write result to file
			with open(ma_datfile+'.json', 'w') as outfile:
				json.dump(result,outfile, indent=4, sort_keys=True, ensure_ascii=True)

			return True
		else:
			sys.exit("not an MS-MS results file - cannot show peptide summary")
	sys.exit("Failure")

def scantitleformat(title):
	#print "detecting :"+title
	returnformat=''
	for rx in config['SCANTITLE']:
		#print rx
		regex= re.compile(config['SCANTITLE'][rx])
		#print regex.pattern
		m=regex.search(title)
		#print m
		if m is not None:
			#print 'hit:'+rx
			returnformat=rx
		
	if not returnformat:
		print "could not detect measurement name extraction format for scan: "+title
		print "aborting..."
		sys.exit(-1)
	
	return returnformat
		
def scantitle(title,scanformat):
	#print title
	regex= re.compile(config['SCANTITLE'][scanformat])
	m=regex.search(title)
	if m is None:
		print "could not extract measurement name with format "+scanformat+" for scan: "+title
		print "aborting..."
		sys.exit(1)
 
	fullm=m.group(1).replace('\\','\/') #replace \ with / (linux path style)
	partm=fullm.split('\/') #split
	
	return partm[-1] #return last element	

def make_modpeptide(ma_mods,peptide,modification,removemods):
	peptide_modified=''
	modificationtypes=[] #contains all modification types
	modifications=[] #contains modification and position in peptide
	misscleaved=0
	nterm=int(modification[0:1])
	if nterm > 0:
		if not removemods:
			#add modification 
			peptide_modified+=ma_mods[nterm]['.'][0] # (.) = N-term modification, only expect one
			modifications.append(['0',ma_mods[nterm]['.'][0],ma_mods[nterm]['.'][1]])
		
		if not ma_mods[nterm]['.'][1]=='None' and not ma_mods[nterm]['.'][1] in modificationtypes:
			#add modification type if necessary
			modificationtypes.append(ma_mods[nterm]['.'][1])
	
	modpos=1
	for aa in peptide:
		setmod = int(modification[modpos:modpos+1])
		if setmod > 0:
			#loop modifications like STY
			for k in ma_mods[setmod]:
				#if amino acid matches key in modification dict
				if k==aa:
					if not removemods:
						#add modification
						peptide_modified+=ma_mods[setmod][k][0]
						modifications.append([modpos,ma_mods[setmod][k][0],ma_mods[setmod][k][1]])
					else:
						#add just the amino acid without modification
						peptide_modified+=aa
					if not ma_mods[setmod][k][1]=='None' and not ma_mods[setmod][k][1] in modificationtypes:
						#add modification type
						modificationtypes.append(ma_mods[setmod][k][1])
					#leave modification loop
					break
		else:
			peptide_modified+=aa
		#check misscleavages
		if modpos < len(peptide):
			#print aa,peptide[modpos:modpos+1] 
			if (aa == "R" or aa == "K") and not peptide[modpos:modpos+1] == "P":
				misscleaved+=1
		#increase position
		modpos+=1
		
	#return modified peptide, modification type, modifications and misscleavages
	return peptide_modified,modificationtypes,modifications,misscleaved
	
if __name__ == '__main__':
	sys.exit(main())

