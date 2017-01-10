import os
from argparse import ArgumentParser
import re
from recup_info import recup_info_star
import subprocess



class Run_Pipeline():
	def __init__(self, innput, output, index):
		self.INPUTPATH=innput
		self.OUTPUTPATH=output
		self.index=index
		#index disponibles: GrCh37 GrCh38
		if index=="GrCh37":
			self.GenecodeANNOTATION="gencode.v25lift37.annotation.gtf"
		elif index=="GrCh38": #pour GrCh38
			self.GenecodeANNOTATION="gencode.v25.annotation.gtf"
		elif index=="Mouse":
			self.GenecodeANNOTATION="gencode.vM11.annotation.gtf"
		else:
			print("error index, choose between Grch37, GrCh38 and Mouse")
			#break
		#DATAPATH="/home/NEXTSEQ/clean_data"

		self.listesample=[]
		for files in os.listdir(self.INPUTPATH):
			if re.match(r"(.)*_R1.fastq.gz$", files):
				element=files.replace("_R1.fastq.gz","")
				self.listesample.append(element)

		#print listesample
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
			try:#os.system('cd '+output)
				os.mkdir(output+'/'+el)
			except:
				pass


	def kmer(self):
		retcode = subprocess.call(["Rscript","Analyse_kmer.R", self.INPUTPATH, self.OUTPUTPATH])

	def STAR(self):
		#RNA-seq aligner
		SCRIPTPATH="~/Documents/RNA-SEQ/STAR"
		STARPATH="/scratch/bin/RNA_Seq/Star"
		if self.index=="GrCh37":
			INDEXSTAR="/scratch/db/genomes/human/index/"
		elif self.index=="GrCh38": #pour GrCh38
			INDEXSTAR="/scratch/db/genomes/human/index/"
		elif self.index=="Mouse":
			#INDEXSTAR="/scratch/db/genomes/mouse/index/"
			INDEXSTAR="/scratch/db/genomes/mus_musculus/index"
		else:
			print("error index, choose between Grch37, GrCh38 and Mouse")
			#break
		
		#TODO changer de place les fichiers d'annot GTF
		for element in self.listesample:
			os.system("/"+STARPATH+"/STAR --genomeDir "+INDEXSTAR+" --readFilesIn "+self.INPUTPATH+"/"+element+"*R2*.fastq.gz "+self.INPUTPATH+"/"+element+"*R1*.fastq.gz --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 18 --chimScoreMin 12 --outReadsUnmapped Fastx --outFileNamePrefix "+self.OUTPUTPATH+"/"+element+"/"+element+" --sjdbGTFfile ~/db/"+self.GenecodeANNOTATION+" --outSAMtype BAM SortedByCoordinate")



	def recup_info_STAR(self):
		####exploitation des resultats de star####
		#regroupement des informations du mapping de chaque echantillon
		liste_info=[]
		for element in self.listesample:
			liste_info.append(self.OUTPUTPATH+"/"+element+"/"+str(element)+"Log.final.out")
		overallstatsFile=str(self.OUTPUTPATH+"/overallstats.txt")
		recup_info_star(liste_info,overallstatsFile)





	def featurecount(self):		
		#######featurecounts  read summarization program
		#os.system("cd ~/Documents/RNA-SEQ/STAR/subread-1.5.1-source/bin/")
		for element in self.listesample:
			os.system("~/Documents/RNA-SEQ/STAR/subread-1.5.1-source/bin/featureCounts -g gene_name -p -s 1 -T 4 -a ~/db/"+self.GenecodeANNOTATION+" -o "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt "+self.OUTPUTPATH+"/"+element+"/"+element+"Aligned.sortedByCoord.out.bam")
		#os.system("cd -")
		


	def group_data(self):
		for element in self.listesample:
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_count.txt; do namefile=`basename $files .txt`; sed 1,2d $files -i; cut -f1,7 $files | sort -k 1b,1 > "+self.OUTPUTPATH+"/"+element+"/$namefile'sorted.txt'; done")

			#os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_countsorted.txt; do sed 1d $files -i; done")


		os.system("awk '{print $1}' "+self.OUTPUTPATH+"/"+self.listesample[0]+"/"+self.listesample[0]+"_countsorted.txt>"+self.OUTPUTPATH+"/allcounts.txt")
		liste_samples_for_allcount="GeneId"
		
		for element in self.listesample:
			liste_samples_for_allcount=liste_samples_for_allcount+" "+element
			os.system("for files in "+self.OUTPUTPATH+"/"+element+"/"+element+"_countsorted.txt; do namefile=`basename $files _countsorted.txt`; sort -k 1b,1 "+self.OUTPUTPATH+"/allcounts.txt>allcountsorted.txt; echo $files; join allcountsorted.txt $files > "+self.OUTPUTPATH+"/allcounts.txt; rm allcountsorted.txt; done")

		os.system("echo "+liste_samples_for_allcount+" | cat - "+self.OUTPUTPATH+"/allcounts.txt > temp && mv temp "+self.OUTPUTPATH+"/allcounts.txt")



if __name__=='__main__':
	
	description = ("input: fastq,  Ouput= count table for differential expression")
	parser = ArgumentParser(description=description)
	parser.add_argument('i', help="entry file with all fastq")
	parser.add_argument('o', help="output directory")
	parser.add_argument('-id','--index', choices=['GrCh37','GrCh38','Mouse'])
	args=parser.parse_args()

	MyRun=Run_Pipeline(args.i, args.o, args.index)
	MyRun.kmer()
	#MyRun.STAR()
	#MyRun.recup_info_STAR()
	#MyRun.featurecount()
	#MyRun.group_data()





