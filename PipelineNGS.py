#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Pipeline de RNA-seq permettant:
	-une pré-analyse de fichiers fastq par analyses des kmer
	-Un mapping sur le génome choisi, via le logiciel STAR
	-groupement par gene-id des reads obtenus lors du mapping, via le logiciel featurecount
	-Analyses d'expression differentielle, par le logiciel DESeq2
"""	


import os
from argparse import ArgumentParser
import re
from recup_info import recup_info_star
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
			elif "FEATURECOUNTSPATH=" in line:
				self.FEATURECOUNTSPATH=line.split("=")[1][:-1]
			elif "SALMONPATH=" in line:
				self.SALMONPATH=line.split("=")[1][:-1]
		
		self.INPUTPATH=innput
		self.OUTPUTPATH=output
		Genecode=False
		indexStar=False
		salmonindex=False
		Trouve1=False
		Trouve2=False
		fichier_config = open(str(fileconfig))
		contenu_config = fichier_config.readlines()
		for line in contenu_config:
			if "####Index STAR####" in line:
				indexStar=True
			if indexStar==True:
				if self.index in line:
					self.INDEXSTAR=line.split("=")[1][:-1]
					indexStar=False
					Trouve1=True
			if "####GENCODE#####" in line:
				Genecode=True
				indexStar=False
			if Genecode==True:
				if self.index in line:
					self.GenecodeANNOTATION=line.split("=")[1][:-1]
					Genecode=False
					Trouve2=True
			if "##INDEX SALMON###" in line:
				salmonindex=True
			if salmonindex==True:
				if self.index in line:
					self.SALMONINDEX=line.split("=")[1][:-1]
					salmonindex=False
		if Trouve1==False or Trouve2==False:
			print("\n\ngenome \""+str(self.index)+"\" non valide\n\n")






		if self.reads=="paired":
			self.listesample=[]
			for files in os.listdir(self.INPUTPATH):
				if re.match(r"(.)*_R1.fastq.gz$", files):
					element=files.replace("_R1.fastq.gz","")
					self.listesample.append(element)
		else:
			self.listesample=[]
			for files in os.listdir(self.INPUTPATH):
				if re.match(r"(.)*.fq.gz$", files):
					element=files.replace(".fq.gz","")
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
		if steps=="all" or steps=="diffexpress":

			try:
				os.mkdir(output+'/DiffExpress')
			except:
				pass
			
			mon_fichier = open(str(fileconfig))
			contenu = mon_fichier.readlines()
			fichier_serie = open(output+"/DiffExpress/serie.txt", "w")
			fichier_test= open(output+"/DiffExpress/compare.txt", "w")
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
						fichier_test= open(output+"/DiffExpress/compare.txt", "a")
				elif Part2==False and Part1==True:
					fichier_serie.write(line)					
				elif Part2==True:
					fichier_test.write(line)
			

	def kmer(self):
		"""
		Analyse des échantillons par kmer (avant mapping), creation d'un premier arbre des echantillons
		"""
		retcode = subprocess.call(["Rscript","Analyse_kmer.R", self.INPUTPATH, self.OUTPUTPATH])

	def STAR(self):
		"""
		Mapping des fastas sur le genome
		logiciel STAR
		"""
		# si paired-end
		if self.reads=="paired":
			for element in self.listesample:
				os.system("/"+self.STARPATH+"/STAR --genomeDir "+self.INDEXSTAR+" --readFilesIn "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 18 --chimScoreMin 12 --outReadsUnmapped Fastx --outFileNamePrefix "+self.OUTPUTPATH+"/"+element+"/"+element+" --sjdbGTFfile "+self.GenecodeANNOTATION+" --outSAMtype BAM SortedByCoordinate")
		# si single-end
		else:
			for element in self.listesample:
				os.system("/"+self.STARPATH+"/STAR --genomeDir "+self.INDEXSTAR+" --readFilesIn "+self.INPUTPATH+"/"+element+".fq.gz --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 18 --chimScoreMin 12 --outReadsUnmapped Fastx --outFileNamePrefix "+self.OUTPUTPATH+"/"+element+"/"+element+" --sjdbGTFfile "+self.GenecodeANNOTATION+" --outSAMtype BAM SortedByCoordinate")


	def recup_info_STAR(self):
		"""
		creation d'un fichier recapitulant les resultats obtenu lors de l'alignement (couverture, reads mappes etc) pour chaque echantillon
		"""
		liste_info=[]
		for element in self.listesample:
			liste_info.append(self.OUTPUTPATH+"/"+element+"/"+str(element)+"Log.final.out")
		overallstatsFile=str(self.OUTPUTPATH+"/overallstats.txt")
		recup_info_star(liste_info,overallstatsFile)





	def featurecount(self):		
		"""
		read summarization program
		utlisation du logiciel featurecounts pour regrouper les resultats du mapping sous forme de tableau de comptage des reads obtenus, par gene_id (possibilite de changer)
		"""
		#Si paired-end
		if self.reads=="paired":
			for element in self.listesample:
				os.system(self.FEATURECOUNTSPATH+"/featureCounts -g gene_id -p -s 1 -T 4 -a "+self.GenecodeANNOTATION+" -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt "+self.OUTPUTPATH+"/"+element+"/"+element+"Aligned.sortedByCoord.out.bam")

		else:
		#Si single-end
			for element in self.listesample:
				os.system(self.FEATURECOUNTSPATH+"/featureCounts -T 4 -t exon -g gene_id -a "+self.GenecodeANNOTATION+" -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt "+self.OUTPUTPATH+"/"+element+"/"+element+"Aligned.sortedByCoord.out.bam")



	def group_data(self):
		"""
		Tri et regroupement  en un seul fichier des resultats obtenus apres featurecounts pour chaque echantillon, qui sera utilisable pour l analyse de l expression differentielle
		"""

		for element in self.listesample:
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt; do namefile=`basename $files .txt`; sed 1,2d $files -i; cut -f1,7 $files | sort -k 1b,1 > "+self.OUTPUTPATH+"/"+element+"/$namefile'sorted.txt'; done")


		os.system("awk '{print $1}' "+self.OUTPUTPATH+"/"+self.listesample[0]+"/"+self.listesample[0]+"_countsorted.txt>"+self.OUTPUTPATH+"/allcounts.txt")
		liste_samples_for_allcount="GeneId"
		
		for element in self.listesample:
			liste_samples_for_allcount=liste_samples_for_allcount+" "+element
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_countsorted.txt; do namefile=`basename $files _countsorted.txt`; sort -k 1b,1 "+self.OUTPUTPATH+"/allcounts.txt>allcountsorted.txt; echo $files; join allcountsorted.txt $files > "+self.OUTPUTPATH+"/allcounts.txt; rm allcountsorted.txt; done")

		os.system("echo "+liste_samples_for_allcount+" | cat - "+self.OUTPUTPATH+"/allcounts.txt > temp && mv temp "+self.OUTPUTPATH+"/allcounts.txt")



	def add_gene_name(self):
		mon_fichier = open(self.OUTPUTPATH+"/allcounts.txt", "r")
		contenu = mon_fichier.readlines()

		mon_fichier2 = open(self.GenecodeANNOTATION[1:-1], "r")
		contenu2 = mon_fichier2.readlines()
		fichier_output= open (self.OUTPUTPATH+"/allcounts_gene.txt", "w")
		fichier_output.write("gene_name"+" "+contenu[0])
		dic_liste_id={}
		for lignes in contenu[1:]:
			lignessplit=lignes.split()
			nameTranscrit=lignessplit[0]
			dic_liste_id[nameTranscrit]=lignes


		for ligne in contenu2:
				pos0= ligne.find('gene_id')
				gene_id=ligne[pos0+8:-1]
				geneid=gene_id.split(";")[0][1:-1]	
				if geneid in dic_liste_id.keys():
		
					pos1 = ligne.find('gene_name')
					gene_name=ligne[pos1+10:-1]
					gene=gene_name.split(";")

					fichier_output.write(gene[0][1:-1]+" "+dic_liste_id.get(geneid))
					del dic_liste_id[geneid]

		fichier_output.close()



		


	def differential_expression(self):
		"""
		analyse de l expression differentielle entre differentes conditions
		fichier de configuration necessaire (voir exemple)
		"""

		subprocess.call(["Rscript","differential_expression.R", self.OUTPUTPATH+"/DiffExpress/", self.OUTPUTPATH+"/allcounts_gene.txt", self.index])
		




	def salmon(self):
		for element in self.listesample:
			if self.reads=="paired":
				os.system(self.SALMONPATH+"/salmon quant -i "+self.SALMONINDEX+" -l A -1 "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz  -2 "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant")
			#else:
				#TODO !!!!!!!!!!!!!!!!!!!!!

		for element in self.listesample:
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/quant.sf; do sed 1,2d $files -i; cut -f1,4 $files | sort -k 1b,1 > "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/"+element+"sorted.txt; done")


		os.system("awk '{print $1}' "+self.OUTPUTPATH+"/"+self.listesample[0]+"/"+self.listesample[0]+"_transcripts_quant/"+self.listesample[0]+"sorted.txt>"+self.OUTPUTPATH+"/all_estimate_counts.txt")
		liste_samples_for_all_estimate_counts="TrancritId"
		
		for element in self.listesample:
			liste_samples_for_all_estimate_counts=liste_samples_for_all_estimate_counts+" "+element
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_transcripts_quant/"+element+"sorted.txt; do namefile=`basename $files sorted.txt`; sort -k 1b,1 "+self.OUTPUTPATH+"/all_estimate_counts.txt>"+self.OUTPUTPATH+"/all_estimate_countssorted.txt; echo $files; join all_estimate_countssorted.txt $files > "+self.OUTPUTPATH+"/all_estimate_counts.txt; rm "+self.OUTPUTPATH+"/all_estimate_countssorted.txt; done")

		os.system("echo "+liste_samples_for_all_estimate_counts+" | cat - "+self.OUTPUTPATH+"/all_estimate_counts.txt > temp && mv temp "+self.OUTPUTPATH+"/all_estimate_counts.txt")






if __name__=='__main__':
	
	description = ("input: fastq,  Ouput= count table for differential expression")
	parser = ArgumentParser(description=description)
	parser.add_argument('-i', help="entry file with all fastq")
	parser.add_argument('-o', help="output directory")
	parser.add_argument('-s','--steps', choices=['all','star','featurecount','diffexpress', 'salmon'])
	parser.add_argument('-c','--config', help="configFile")
	args=parser.parse_args()

	MyRun=Run_Pipeline(args.i, args.o, args.steps, args.config)
	if args.steps=="all" or args.steps=="star":
		MyRun.STAR()
		MyRun.recup_info_STAR()
	if args.steps=="all" or args.steps=="featurecount":
		#MyRun.featurecount()
		#MyRun.group_data()
		MyRun.add_gene_name()
	if args.steps=="all" or args.steps=="diffexpress":
		MyRun.differential_expression()
	if args.steps=="salmon":
		MyRun.salmon()

