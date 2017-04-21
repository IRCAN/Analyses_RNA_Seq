import sys


def recup_info_star(listesample,overallstatsFile):
	liste_sample = listesample
	allsample=""
	#listes
	nbr_input_reads="Number_of_input_reads"
	reads_length="Average_input_read_length"
	uniquely_mapped_reads="%_Uniquely_mapped_reads"
	reads_mapped_multiple_loci="%_of_reads_mapped_to_multiple_loci"
	#mapped_too_many_loci=""
	unmapped_bc_mismatches="%_of_reads_unmapped:_too_many_mismatches"
	unmapped_bc_short="%_of_reads_unmapped:_too_short"
	unmapped_bc_other="%_of_reads_unmapped:_other"


	for sample in liste_sample:
		allsample=allsample+"\t"+sample.split("/")[-2]
		with open(sample, "r") as file_sample:
			contenu = file_sample.readlines()
			for lignes in contenu:
				if "Number of input reads" in lignes:
					texte=lignes.split("\t")
					nbr_input_reads=nbr_input_reads+"\t"+(texte[-1][:-1])
				elif "Average input read length" in lignes:
					texte=lignes.split("\t")
					reads_length=reads_length+"\t"+(texte[-1][:-1])
				elif "Uniquely mapped reads %" in lignes:
					texte=lignes.split("\t")
					uniquely_mapped_reads=uniquely_mapped_reads+"\t"+(texte[-1][:-2]) #[:-2] pour retirer "%"
				elif "% of reads mapped to multiple loci" in lignes:
					texte=lignes.split("\t")
					reads_mapped_multiple_loci=reads_mapped_multiple_loci+"\t"+(texte[-1][:-2])
				#elif "% of reads mapped to too many loci" in lignes:
				#	texte=lignes.split("\t")
				#	mapped_too_many_loci.append(texte[-1][:-1])
				elif "% of reads unmapped: too many mismatches" in lignes:
					texte=lignes.split("\t")
					unmapped_bc_mismatches=unmapped_bc_mismatches+"\t"+(texte[-1][:-2])
				elif "% of reads unmapped: too short" in lignes:
					texte=lignes.split("\t")
					unmapped_bc_short=unmapped_bc_short+"\t"+(texte[-1][:-2])
				elif " % of reads unmapped: other" in lignes:
					texte=lignes.split("\t")
					unmapped_bc_other=unmapped_bc_other+"\t"+(texte[-1][:-2])


	with open(overallstatsFile, "w") as overall_stats:
		overall_stats.write(allsample+"\n")
		overall_stats.write(nbr_input_reads+"\n")
		overall_stats.write(reads_length+"\n")
		overall_stats.write(uniquely_mapped_reads+"\n")
		overall_stats.write(reads_mapped_multiple_loci+"\n")
		#overall_stats.write(str(mapped_too_many_loci)[1:-1]+"\n")
		overall_stats.write(unmapped_bc_mismatches+"\n")
		overall_stats.write(unmapped_bc_short+"\n")
		overall_stats.write(unmapped_bc_other+"\n")
