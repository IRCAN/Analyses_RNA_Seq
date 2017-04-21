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
	
	def __init__(self, innput, output, index, reads, *DiffExpress):
		"""
		Recupere les fichiers fasta a analyser
		Recupere le fichier d annotation requis
		Cree les repertoires de sorties si non existants
	
		"""
		self.INPUTPATH=innput
		self.OUTPUTPATH=output
		self.index=index
		self.reads=reads
		#index disponibles:  TODO: changer repertoire des fichiers annotes
		if index=="GrCh37":
			self.GenecodeANNOTATION="/home/Share/ftessier/PipelineNGS/gencode.v25lift37.annotation.gtf"
		elif index=="GrCh38": 
			self.GenecodeANNOTATION="/home/Share/ftessier/PipelineNGS/gencode.v25.annotation.gtf"
		elif index=="Mouse":
			self.GenecodeANNOTATION="/home/Share/ftessier/PipelineNGS/gencode.vM11.annotation.gtf"
		elif index=="Mm9":
			self.GenecodeANNOTATION="/home/Share/ftessier/PipelineNGS/gencode.vM1.annotation.gtf" 
		else:
			print("error index, choose between Grch37, GrCh38 and Mouse")


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
		if DiffExpress[1] or DiffExpress[0] is not None:

			try:
				os.mkdir(output+'/DiffExpress')
			except:
				pass
			
		        Part2=False
			try:
				mon_fichier = open(str(DiffExpress[1]))
			except:
				mon_fichier = open(str(DiffExpress[0]))
			contenu = mon_fichier.readlines()
			fichier_serie = open(output+"/DiffExpress/serie.txt", "w")
			fichier_test= open(output+"/DiffExpress/compare.txt", "w")		
			for line in contenu:
				if line[0:2] != "--" and Part2==False:
					fichier_serie.write(line)
				elif Part2==False:
					Part2=True
					fichier_serie.close()
				else:
					fichier_test= open(output+"/DiffExpress/compare.txt", "a")
					fichier_test.write(line)
			fichier_test.close()

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

		STARPATH="/scratch/bin/RNA_Seq/Star"
		if self.index=="GrCh37":
			INDEXSTAR="/scratch/db/genomes/human/index/GrCh37_chr_only"
		elif self.index=="GrCh38": 
			INDEXSTAR="/scratch/db/genomes/human/index/GrCh38_chr_only"
		elif self.index=="Mouse":
			INDEXSTAR="/scratch/db/genomes/mus_musculus/index/mm10"
		elif self.index=="Mm9":

			INDEXSTAR="/scratch/db/genomes/mus_musculus/index/mm9"
		else:
			print("error index, choose between Grch37, GrCh38 and Mouse")

		# si paired-end
		if self.reads=="paired":
			for element in self.listesample:
				os.system("/"+STARPATH+"/STAR --genomeDir "+INDEXSTAR+" --readFilesIn "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 18 --chimScoreMin 12 --outReadsUnmapped Fastx --outFileNamePrefix "+self.OUTPUTPATH+"/"+element+"/"+element+" --sjdbGTFfile "+self.GenecodeANNOTATION+" --outSAMtype BAM SortedByCoordinate")
		# si single-end
		else:
			for element in self.listesample:
				os.system("/"+STARPATH+"/STAR --genomeDir "+INDEXSTAR+" --readFilesIn "+self.INPUTPATH+"/"+element+".fq.gz --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 18 --chimScoreMin 12 --outReadsUnmapped Fastx --outFileNamePrefix "+self.OUTPUTPATH+"/"+element+"/"+element+" --sjdbGTFfile "+self.GenecodeANNOTATION+" --outSAMtype BAM SortedByCoordinate")


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
				os.system("/scratch/bin/alignment/subread/bin/featureCounts -g gene_id -p -s 1 -T 4 -a "+self.GenecodeANNOTATION+" -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt "+self.OUTPUTPATH+"/"+element+"/"+element+"Aligned.sortedByCoord.out.bam")

		else:
		#Si single-end
			print(self.listesample)
			for element in self.listesample:
				os.system("/scratch/bin/alignment/subread/bin/featureCounts -T 4 -t exon -g gene_id -a "+self.GenecodeANNOTATION+" -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt "+self.OUTPUTPATH+"/"+element+"/"+element+"Aligned.sortedByCoord.out.bam")



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

		mon_fichier2 = open(self.GenecodeANNOTATION, "r")
		contenu2 = mon_fichier2.readlines()
		fichier_output= open ("allcounts_gene.txt", "w")
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



		


	def differential_expression(self, index):
		"""
		analyse de l expression differentielle entre differentes conditions
		fichier de configuration necessaire (voir exemple)
		"""

		subprocess.call(["Rscript","differential_expression.R", self.OUTPUTPATH+"/DiffExpress/", self.OUTPUTPATH+"/allcounts_gene.txt", index])





if __name__=='__main__':
	
	description = ("input: fastq,  Ouput= count table for differential expression")
	parser = ArgumentParser(description=description)
	parser.add_argument('i', help="entry file with all fastq")
	parser.add_argument('o', help="output directory")
	parser.add_argument('-id','--index', choices=['GrCh37','GrCh38','Mouse','Mm9'])
	parser.add_argument('-seq','--sequencing', choices=['paired','single'])
	parser.add_argument("--Kmer")
	parser.add_argument("--Star")
	parser.add_argument("--Featurecount")
	parser.add_argument("--DiffExpress")
	parser.add_argument("--All", help="Run all the options since the beginning")
	args=parser.parse_args()

	MyRun=Run_Pipeline(args.i, args.o, args.index, args.sequencing, args.All, args.DiffExpress)
	if args.Kmer is not None or args.All is not None:
		MyRun.kmer()
	if args.Star is not None or args.All is not None:
		MyRun.STAR()
		MyRun.recup_info_STAR()
	if args.Featurecount is not None or args.All is not None:
		MyRun.featurecount()
		MyRun.group_data()
		MyRun.add_gene_name()
	if args.DiffExpress is not None or args.All is not None:
		MyRun.differential_expression( args.index)



