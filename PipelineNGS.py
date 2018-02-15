#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Pipeline de RNA-seq permettant:
	-une pré-analyse de fichiers fastq par analyses des kmer
	-Un mapping sur le génome choisi, via le logiciel STAR
	-ou un pseudo-alignement sur le transcriptome, avec le logiciel Salmon
	-groupement par gene-id des reads obtenus lors du mapping, via le logiciel featurecount
	-Analyses d'expression differentielle, par le logiciel DESeq2
"""	


import os
from argparse import ArgumentParser
import re
import subprocess



class Run_Pipeline():

	def __init__(self, innput, output, steps, fileconfig):
		"""
		Recupere les fichiers fasta a analyser
		Recupere le fichier d annotation requis
		Cree les repertoires de sorties si non existants
	
		"""
		fichier_config = open(str(fileconfig))
		contenu_config = fichier_config.readlines()
		for line in contenu_config:
			if ">Paired-end" in line:
				if "TRUE" in line:
					self.reads="paired"
				else:
					self.reads="single"
			elif ">Genome" in line:
				self.index=line.split("=")[1][:-1]
			elif "STARPATH=" in line:
				self.STARPATH=line.split("=")[1][:-1]
			elif "SALMONPATH=" in line:
				self.SALMONPATH=line.split("=")[1][:-1]
			elif "RSEMPATH=" in line:
				self.RSEMPATH=line.split("=")[1][:-1]
		self.INPUTPATH=innput
		self.OUTPUTPATH=output
		Genecode=False
		indexStar=False
		salmonindex=False
		kallistoindex=False
		rsemindex=False
		tx2gene=False

		fichier_config = open(str(fileconfig))
		contenu_config = fichier_config.readlines()
		for line in contenu_config:
			if "####GENCODE#####" in line:
				Genecode=True
				indexStar=False
			if Genecode==True:
				if self.index in line:
					self.GenecodeANNOTATION=line.split("=")[1][:-1]
					Genecode=False
			if "##INDEX SALMON###" in line:
				salmonindex=True
			if salmonindex==True:
				if self.index in line:
					self.SALMONINDEX=line.split("=")[1][:-1]
					salmonindex=False
			if "####Index RSEM##" in line:
				rsemindex=True
			if rsemindex==True:
				if self.index in line:
					self.RSEMINDEX=line.split("=")[1][:-1]
					rsemindex=False
			if "####TX2GENE###" in line:
				tx2gene=True
			if tx2gene==True:
				if self.index in line:
					self.TX2GENE=line.split("=")[1][:-1]
					tx2gene=False







		if self.reads=="paired":
			self.listesample=[]
			for files in os.listdir(self.INPUTPATH):
				if re.match(r"(.)*_R1.fastq.gz$", files):
					element=files.replace("_R1.fastq.gz","")
					self.listesample.append(element)
				elif re.match(r"(.)*_R1.fastq", files):
					element=files.replace("_R1.fastq","")
					self.listesample.append(element)
		else:
			self.listesample=[]
			for files in os.listdir(self.INPUTPATH):
				if re.match(r"(.)*.fastq", files):
					element=files.replace(".fastq","")
					self.listesample.append(element)
		try:
			self.listesample.remove("Undetermined_S0")
		except:
			pass
		#if necessary, make the output directory
		try:
			os.mkdir(output)
		except:
			pass
		for el in self.listesample:
			try:
				os.mkdir(output+'/'+el)
			except:
				pass

		if steps=="salmon":

			try:
				os.mkdir(output+'/DiffExpressKmer')
			except:
				pass
			
			mon_fichier = open(str(fileconfig))
			contenu = mon_fichier.readlines()
			fichier_serie = open(output+"/DiffExpressKmer/serie.txt", "w")
			fichier_test= open(output+"/DiffExpressKmer/compare.txt", "w")
			Part1=False
			Part2=False		
			for line in contenu:
				if line[0:2] == "--":
					if Part1==False and Part2==False:
						Part1=True
					elif Part2==True:
						Part2=False
						fichier_test.close()
					else:
						Part2=True
						Part1=False
						fichier_serie.close()
						fichier_test= open(output+"/DiffExpressKmer/compare.txt", "a")
				elif Part2==False and Part1==True:
					fichier_serie.write(line)					
				elif Part2==True:
					fichier_test.write(line)
		


		if steps=="rsem" :

			try:
				os.mkdir(output+'/DiffExpressRSEM')
			except:
				pass
			
			mon_fichier = open(str(fileconfig))
			contenu = mon_fichier.readlines()
			fichier_serie = open(output+"/DiffExpressRSEM/serie.txt", "w")
			fichier_test= open(output+"/DiffExpressRSEM/compare.txt", "w")
			Part1=False
			Part2=False		
			for line in contenu:
				if line[0:2] == "--":
					if Part1==False and Part2==False:
						Part1=True
					elif Part2==True:
						Part2=False
						fichier_test.close()
					else:
						Part2=True
						Part1=False
						fichier_serie.close()
						fichier_test= open(output+"/DiffExpressRSEM/compare.txt", "a")
				elif Part2==False and Part1==True:
					fichier_serie.write(line)					
				elif Part2==True:
					fichier_test.write(line)


		#END############preparation of files and directories according to the selected options	##################################################################
###################################################################################################################################################################################################

	def kmer(self):
		"""
		Analyse des échantillons par kmer (avant mapping), creation d'un premier arbre des echantillons
		"""
		retcode = subprocess.call(["Rscript","/home/Share/ftessier/PipelineNGS/Analyse_kmer.R", self.INPUTPATH, self.OUTPUTPATH])

###################################################################################################################################################################################################


	def rsem(self):
		"""
		mapping par star, via RSEM
		"""
		for element in self.listesample:
				os.system("/"+self.RSEMPATH+"/rsem-calculate-expression --star --star-path "+self.STARPATH+" --star-gzipped-read-file -p 4 --paired-end "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz "+self.RSEMINDEX+" "+self.OUTPUTPATH+"/"+element+"/"+element)


###################################################################################################################################################################################################

	def fix_Id_for_gene_name(self):
		for element in self.listesample:
			try:
				os.system("cp "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.old.sf "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.sf")
			except:
				pass
			try:
				os.system("cp "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.sf "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.old.sf")
			except:
				pass
			#mon_fichier = open(self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/abundance.old.tsv", "r")

			mon_fichier = open(self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.old.sf", "r")

			mon_fichier_output = open(self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.sf", "w")
			contenu = mon_fichier.readlines()

			for lignes in contenu[0]:
				mon_fichier_output.write(lignes)
			for lignes in contenu[1:]:
				l=lignes.split("	", 1)
				a=l[0].split(".")
				mon_fichier_output.write(str(a[0])+"	"+str("	".join(l[-1:])))

###################################################################################################################################################################################################

	def salmon(self):
		for element in self.listesample:
			if self.reads=="paired":
				os.system(self.SALMONPATH+"/salmon quant -i "+self.SALMONINDEX+" -l A -1 "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz  -2 "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz --numBootstraps 100 -p 8 -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant -g "+self.GenecodeANNOTATION)
			else:
				os.system(self.SALMONPATH+"/salmon quant -i "+self.SALMONINDEX+" -l A -r "+self.INPUTPATH+"/"+element+".fastq --numBootstraps 100 -p 4 -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant")
		
###################################################################################################################################################################################################

	def differential_gene_expression(self, step):
		"""
		analyse de l expression differentielle entre differentes conditions
		fichier de configuration necessaire (voir exemple)
		"""
		subprocess.call(["Rscript","/home/Share/ftessier/PipelineNGS/differential_gene_expression.R", self.OUTPUTPATH, self.index, step , self.TX2GENE])

###################################################################################################################################################################################################

	def differential_transcript_expression(self, step):
		subprocess.call(["Rscript","/home/Share/ftessier/PipelineNGS/Differential_transcript_expression.R", self.OUTPUTPATH, self.index, step, self.TX2GENE])


###################################################################################################################################################################################################	
############à modifier############
	def ajout_legend(self): 
		os.system("for file in "+self.OUTPUTPATH+"/DiffExpress/*/*_*regulate*.png; do repertoire=$(dirname $file); condition=$(basename $repertoire); python /home/Share/ftessier/PipelineNGS/legend_image.py $file $repertoire/ALL_$condition.xls "+self.OUTPUTPATH+"/DiffExpress/$condition/keggresgreater.txt "+self.OUTPUTPATH+"/DiffExpress/$condition; done ; mkdir "+self.OUTPUTPATH+"/DiffExpress/$condition/Pathways_Upregulate ; mv "+self.OUTPUTPATH+"/DiffExpress/$condition/Qvalue* "+self.OUTPUTPATH+"/DiffExpress/$condition/Pathways_Upregulate")

		os.system("for file in "+self.OUTPUTPATH+"/DiffExpress/*/*_*downregulate*.png; do repertoire=$(dirname $file); condition=$(basename $repertoire); python /home/Share/ftessier/PipelineNGS/legend_image.py $file $repertoire/ALL_$condition.xls "+self.OUTPUTPATH+"/DiffExpress/$condition/keggresless.txt "+self.OUTPUTPATH+"/DiffExpress/$condition; done; mkdir "+self.OUTPUTPATH+"/DiffExpress/$condition/Pathways_Downregulate ; mv "+self.OUTPUTPATH+"/DiffExpress/$condition/Qvalue* "+self.OUTPUTPATH+"/DiffExpress/$condition/Pathways_Downregulate")
		os.system("for file in "+self.OUTPUTPATH+"/DiffExpress/*/*.xml; do rm $file; done")
		os.system("for file in "+self.OUTPUTPATH+"/DiffExpress/*/*.png; do rm $file; done")


###################################################################################################################################################################################################




if __name__=='__main__':
	
	description = ("...")
	parser = ArgumentParser(description=description)
	parser.add_argument('-i', help="entry directory with all fastq")
	parser.add_argument('-o', help="output directory")
	parser.add_argument('-s','--steps', choices=['kmer','diffexpress', 'salmon', 'rsem'])
	parser.add_argument('-c','--config', help="configFile")
	args=parser.parse_args()

	MyRun=Run_Pipeline(args.i, args.o, args.steps, args.config)
	if args.steps=="kmer":
		MyRun.kmer()
	if args.steps=="salmon" :
		MyRun.salmon()
		MyRun.fix_Id_for_gene_name()
	if args.steps=="rsem" :
		MyRun.rsem()
	if args.steps=="rsem" or args.steps=="salmon":
		MyRun.differential_gene_expression(args.steps)
		MyRun.differential_transcript_expression(args.steps)
		MyRun.ajout_legend()
