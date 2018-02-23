#!/usr/bin/env python
"""
p2pmysql.py: pep2pro mysql database class

"""
import mysql.connector	
import json
import sys

class pep2prodb:
######################################################################################		
	def __init__(self,user,password,host,database):
		cnx = mysql.connector.connect(user=user
				,password=password
				,host=host
				,database=database)	
		self.cnx=cnx #hold database connection
		self.modificationtypes={} #holds all modificationtypes
		self.modifications={} #holds modifications
		self.proteins={}
		self.peptides={}
		self.peptides_modified={}
		self.measurements={}
######################################################################################		
	def close():
		self.cnx.close()
######################################################################################		
	def mysql_commit(self):
		self.cnx.commit()
######################################################################################		
	def mysql_rollback(self):
		self.cnx.rollback()
######################################################################################		
	def mysql_execute(self,sql,params,commit=0):
		csql = self.cnx.cursor()
		#print sql,params
		try:
			csql.execute(sql,params)
			if commit==1:
				self.cnx.commit()
			return csql
		except mysql.connector.Error as err:
			print("Error: Executing SQL Command: {}".format(err))	
			print "-"*60
			print csql.statement
			print "-"*60
			raise
######################################################################################		
	def get_searchdatabase_id(self,searchdb):
		searchdbid=None
		sql = ("select id from searchdatabase "
				" where searchdb = %s")
		csr=self.mysql_execute(sql,[searchdb])
		rows=csr.fetchall()
		for ([id]) in rows:
			searchdbid=id
		return searchdbid
######################################################################################		
	def create_searchdatabase_id(self,searchdb,searchdatabase):
		sql = ("insert into searchdatabase "
				" (searchdb,searchdatabase) "
				" values(%s,%s)")
		csr=self.mysql_execute(sql,[searchdb,searchdatabase])
		searchdb=csr.lastrowid
		return searchdb
######################################################################################		
	def load_searchdatabase(self,config,searchdbid,searchdatabase):
		with open(config['SEARCHDATABASEPATH']['PATH']+"/"+searchdatabase+".stats.json") as data_file:    
			data = json.load(data_file)	
			sql = "lock tables protein write"
			csr=self.mysql_execute(sql,[])
			for protident in data:
				sql = ( "insert into protein (searchdatabase_id,accession,locus "
						",description,tryptic_peptides,decoy,protein_length "
						",isoelectric_point,molecular_weight ) "
						" values(%s,%s,%s,%s,%s,%s,%s,%s,%s) ")
				record=[searchdbid,protident,data[protident]['locus']
				,data[protident]['description'],data[protident]['ttp']
				,data[protident]['decoy'],data[protident]['len'] 
				,data[protident]['iep'],data[protident]['mol']
				]
				csr=self.mysql_execute(sql,record)
			sql = "unlock tables"
			csr=self.mysql_execute(sql,[])
			return True
######################################################################################		
	def get_proteins(self,searchdbid):
		sql = ("select a.id,a.accession,a.decoy "
			" from protein a "
			" where a.searchdatabase_id = %s")
		csr=self.mysql_execute(sql,[searchdbid])
		rows=csr.fetchall()
		for ([protein_id,accession,decoy]) in rows:
			self.proteins[accession]={'id':protein_id,'decoy':decoy}
######################################################################################		
	def get_datfile_id(self,datfile,ionsscore,probability,searchdatabase_id,experiment,genotype,tissue,treatment,replicate,fraction):
		datfile_id=False
		#lock table
		sql = "lock tables datfile write"
		csr=self.mysql_execute(sql,[])
		#check if datafile exists and abort
		sql = ("select id from datfile "
				" where datfile=%s "
				" and ionsscore=%s " 
				" and probability=%s "
				" and searchdatabase_id=%s " )
		csr=self.mysql_execute(sql,[datfile,ionsscore,probability,searchdatabase_id])
		rows=csr.fetchall()
		for ([id]) in rows:
			print("Error: datfile already imported.")
			print "-"*60
			sys.exit(1)	
		#insert datfile data
		sql = ("insert into datfile (datfile,ionsscore,probability,searchdatabase_id "
				",experiment,genotype,tissue,treatment,replicate,fraction) "
				" values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s) ")
		record=[datfile,ionsscore,probability,searchdatabase_id
		,experiment,genotype,tissue,treatment,replicate,fraction
		]
		csr=self.mysql_execute(sql,record,1)
		datfile_id=csr.lastrowid
		#unlock table
		sql = "unlock tables"
		csr=self.mysql_execute(sql,[])
		#return id
		return datfile_id
######################################################################################		
	def update_modificationtype(self,modificationtype):
		#check if modificationtype exists or create it
		sql = ("select id,modificationtype "
				" from modificationtype "
				" where modificationtype=%s")
		csr=self.mysql_execute(sql,[modificationtype])
		rows=csr.fetchall()
		for ([id,mtype]) in rows:
			self.modificationtypes[mtype]=id
		if not modificationtype in self.modificationtypes: 
			#create modification data
			sql = ("insert into modificationtype (modificationtype) "
				" values(%s) ")
			csr=self.mysql_execute(sql,[modificationtype],1)
			self.modificationtypes[modificationtype]=csr.lastrowid
######################################################################################		
	def update_modification(self,modification,modificationtype):
		#check if modification exists or create it
		sql = ("select id,modification "
				" from modification "
				" where modification=%s")
		csr=self.mysql_execute(sql,[modification])
		rows=csr.fetchall()
		for ([id,mod]) in rows:
			self.modifications[mod]=id
		if not modification in self.modifications: 
			#create modification data
			sql = ("insert into modification (modification,modificationtype_id) "
				" values(%s,%s) ")
			csr=self.mysql_execute(sql,[modification,self.modificationtypes[modificationtype]],1)
			self.modifications[modification]=csr.lastrowid
######################################################################################		
	def update_peptide(self,peptide,missed_cleavages):
		sql = "lock tables peptide write"
		csr=self.mysql_execute(sql,[])
		#check if datafile exists and abort
		sql = ("select id,peptide "
				" from peptide "
				" where peptide=%s")
		csr=self.mysql_execute(sql,[peptide])
		rows=csr.fetchall()
		for ([id,pep]) in rows:
			self.peptides[pep]=id
		if not peptide in self.peptides: 
			#create modification data
			sql = ("insert into peptide (peptide,missed_cleavages) "
				" values(%s,%s) ")
			csr=self.mysql_execute(sql,[peptide,missed_cleavages],1)
			self.peptides[peptide]=csr.lastrowid
		#unlock table
		sql = "unlock tables"
		csr=self.mysql_execute(sql,[])
######################################################################################		
	def update_peptide_modified(self,peptide_modified,peptide,modifications):
		sql = "lock tables peptide_modified write, modificationtype write, modification write, peptide_modification write"
		csr=self.mysql_execute(sql,[])
		#check if datafile exists and abort
		sql = ("select id,peptide_modified "
				" from peptide_modified "
				" where peptide_modified=%s")
		csr=self.mysql_execute(sql,[peptide_modified])
		rows=csr.fetchall()
		for ([id,pep]) in rows:
			self.peptides_modified[pep]=id
		if not peptide_modified in self.peptides_modified: 
			#create modification data
			sql = ("insert into peptide_modified (peptide_modified,peptide_id) "
				" values(%s,%s) ")
			csr=self.mysql_execute(sql,[peptide_modified,self.peptides[peptide]],1)
			self.peptides_modified[peptide_modified]=csr.lastrowid
			#modificationdetails
			#check each modification and create if necessary
			if not modifications == None:
				for position,modification,modificationtype in modifications:
					#for each modification 
					if not modificationtype in self.modificationtypes:
						self.update_modificationtype(modificationtype)
					if not modification in self.modifications:
						self.update_modification(modification,modificationtype)
					#insert modification details
					sql = ("insert into peptide_modification "
							" (peptide_modified_id,modpos,modification_id) "
							" values(%s,%s,%s) ")
					csr=self.mysql_execute(sql,[self.peptides_modified[peptide_modified],position,self.modifications[modification]],1)
		#unlock table
		sql = "unlock tables"
		csr=self.mysql_execute(sql,[])
#####################################################################################		
	def update_accession(self,accession,peptide):
		rc=False
		sql = "lock tables peptide_protein write"
		csr=self.mysql_execute(sql,[])
		#check if datafile exists and abort
		sql = ("select peptide_id,protein_id "
				" from peptide_protein "
				" where peptide_id=%s"
				" and protein_id=%s")
		csr=self.mysql_execute(sql,[self.peptides[peptide],self.proteins[accession]['id']])
		rows=csr.fetchall()
		for ([pep_id,prot_id]) in rows:
			rc=True
		if not rc: 
			#create modification data
			sql = ("insert into peptide_protein "
					" (peptide_id,protein_id) "
					" values(%s,%s) ")
			csr=self.mysql_execute(sql,[self.peptides[peptide],self.proteins[accession]['id']],1)
		#unlock table
		sql = "unlock tables"
		csr=self.mysql_execute(sql,[])
#####################################################################################		
	def update_measurement(self,measurement):
		sql = "lock tables measurement write"
		csr=self.mysql_execute(sql,[])
		#check if datafile exists and abort
		sql = ("select id,measurement "
				" from measurement "
				" where measurement=%s")
		csr=self.mysql_execute(sql,[measurement])
		rows=csr.fetchall()
		for ([id,mes]) in rows:
			self.measurements[mes]=id
		if not measurement in self.measurements: 
			#create modification data
			sql = ("insert into measurement (measurement) "
				" values(%s) ")
			csr=self.mysql_execute(sql,[measurement],1)
			self.measurements[measurement]=csr.lastrowid
		#unlock table
		sql = "unlock tables"
		csr=self.mysql_execute(sql,[])
		return self.measurements[measurement]
######################################################################################		
	def get_pepetide_peptide_modified_id(self,accessions,peptide,misscleaved,peptide_modified,modifications):
		peptide_id=False
		peptide_modified_id=False
		#check and create peptide
		if not peptide in self.peptides:
			self.update_peptide(peptide,misscleaved)
		#check and create peptide_modified
		if not peptide_modified in self.peptides_modified:
			self.update_peptide_modified(peptide_modified,peptide,modifications)
		#check and create peptide-protein link
		for accession in accessions:
			self.update_accession(accession['protident'],peptide)
		#return ids
		return self.peptides[peptide],self.peptides_modified[peptide_modified]
######################################################################################		
	def update_hit_modificationtypes(self,hit_id,hit_modifications):
		for hitmod in hit_modifications:
			if not hitmod in self.modificationtypes:
				self.update_modificationtype(hitmod)
			sql = ("insert into hit_modificationtype "
					" (hit_id,modificationtype_id) "
					" values(%s,%s) ")
			csr=self.mysql_execute(sql,[hit_id,self.modificationtypes[hitmod]],1)
######################################################################################		

	