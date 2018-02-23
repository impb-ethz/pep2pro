#!/usr/bin/env python
"""
makep2pdb.py: pep2pro script to create search database index,
              protein analysis data and creates an in silico 
              tryptic digest
"""
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Development"

import p2pconfig
import sys,os
import json
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
#argumentparser
import click

config=p2pconfig.read_config()

@click.command()
#which arguments and options are needed? 
#arguments=mandatory, option=possible
#*************** ARGUMENETS *****************************
#input-file
@click.argument('inputfile', type=click.Path(exists=True,readable=True))
def main(inputfile):
	#extract basename
	searchdb=os.path.basename(inputfile)
	#check search db based on basename
	dbconfig=p2pconfig.chk_searchdb(config,searchdb,False) #check searchdb but not json
	#get regular expression for search db
	regexes=p2pconfig.extract_regex(config,dbconfig)
	#define proteins dictionary
	proteins={}
	panalysis={}
	for record in SeqIO.parse(inputfile, 'fasta'):
		#define locus and decoy flag	
		locus = p2pconfig.makelocus(regexes,record.id)
		#check if valid locus 
		if locus['locus'] is not None:
			panalysis[record.id]={}
			panalysis[record.id]['locus']=locus['locus']
			panalysis[record.id]['description']=record.description
			panalysis[record.id]['decoy']=locus['decoy']
			panalysis[record.id]['mol']=-1
			panalysis[record.id]['iep']=-1
			panalysis[record.id]['len']=-1
			panalysis[record.id]['ttp']=-1
			#do sequence analysis if not decoy
		 	if locus['decoy']==False:
				#print locus
				if not str(record.seq) in proteins:
					proteins[str(record.seq)]=[record.id]
				else:
					proteins[str(record.seq)].append(record.id)
				#analysis protein
				analysed_seq = ProteinAnalysis(seq_weight_corrections(str(record.seq)))
				panalysis[record.id]['mol']=analysed_seq.molecular_weight()
				panalysis[record.id]['iep']=analysed_seq.isoelectric_point()
				panalysis[record.id]['len']=len(record.seq)
				
				panalysis[record.id]['ttplist']=list(set(tryptic_digest(str(record.seq))))
				panalysis[record.id]['ttp']=len(panalysis[record.id]['ttplist'])
				del panalysis[record.id]['ttplist'] # not needed anymore 
		
	#restructure data
	jsondata={}
	for seq in proteins:
		for p in proteins[seq]:
			jsondata[p]=proteins[seq]
	#dump data
	with open(inputfile+'.json', 'w') as outfile:
		json.dump(jsondata,outfile, indent=4, sort_keys=True, ensure_ascii=True)
	outfile.close()
	with open(inputfile+'.stats.json', 'w') as outfile:
		json.dump(panalysis,outfile, indent=4, sort_keys=True, ensure_ascii=True)
	outfile.close()

def tryptic_digest(sequence):
	ttp=[]
	cut_sites=[0]
    #define cut-sites
	for i in range(0,len(sequence)-1):
		if sequence[i]=='K' and sequence[i+1]!='P':
			cut_sites.append(i+1)
		elif sequence[i]=='R' and sequence[i+1]!='P':
			cut_sites.append(i+1)
    #add last position if not already done
	if cut_sites[-1]!=len(sequence):
		cut_sites.append(len(sequence))

	if len(cut_sites)>2:
		for j in range(0,len(cut_sites)-1):
			peptide=sequence[cut_sites[j]:cut_sites[j+1]]
			if checksize(peptide):
				#print '->'+peptide
				ttp.append(peptide)
				#subdigest?
				ttp.extend(sub_digest(peptide))
	else: #there is no trypsin site in the protein sequence
		if checksize(sequence):
			#print '->'+peptide
			ttp.append(sequence)	
			#subdigest?
			ttp.extend(sub_digest(sequence))
	return ttp

def seq_weight_corrections(sequence,upperbounds=False):
	#this function is to calculate molecular weight and substitutes
	#ambiguous single letter amino acids by its upper or lower bounds
	if upperbounds:
		sequence=sequence.replace('X','W') # unknown to Tryptophan
		sequence=sequence.replace('B','D') # Aspartic Acid, Asparagine to Asparagine
		sequence=sequence.replace('J','I') # Leucine, Isoleucine to Isoleucine
		sequence=sequence.replace('Z','E') # Glutamic Acid, Glutamine to Glutamine
		
	else:
		sequence=sequence.replace('X','G') # unknown to Glycine
		sequence=sequence.replace('B','N') # Aspartic Acid, Asparagine to Aspartic Acid
		sequence=sequence.replace('J','L') # Leucine, Isoleucine to Leucine
		sequence=sequence.replace('Z','Q') # Glutamic Acid, Glutamine to Glutamic
	return sequence

def sub_digest(peptide):
	ttp=[]
	cut_sites=[0]
	for i in range(0,len(peptide)-1):
		if peptide[i]=='K' or peptide[i]=='R':
			cut_sites.append(i+1)
	
	if len(cut_sites)>1:	
		#print '-->'+peptide
		for j in range(0,len(cut_sites)-1):
			peptide1=peptide[0:cut_sites[j+1]]
			if checksize(peptide1):
				ttp.append(peptide1)
			peptide2=peptide[cut_sites[j+1]:]
			if checksize(peptide2):
				ttp.append(peptide2)
			#print '--->'+peptide1,peptide2
	return ttp

def checksize(peptide):
	#only peptides with a minimal length of 6 amino acids 
	#and a molecular weigth between 400 and 6000 are accepted
	if len(peptide)>=6:
		analysed_seq = ProteinAnalysis(seq_weight_corrections(peptide))
		mol_weight = analysed_seq.molecular_weight()
		if mol_weight >=400 and mol_weight <= 6000:
			return True
		else:
			return False
	else:
		return False	
	
if __name__ == '__main__':
  main()




