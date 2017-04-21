
			
mon_fichier = open("C4_vs_C3.xls", "r")
contenu = mon_fichier.readlines()
mon_fichier2 = open("/home/Share/ftessier/PipelineNGS/gencode.vM11.annotation.gtf", "r")
contenu2 = mon_fichier2.readlines()
fichier_output= open ("C4vsC3.txt", "w")
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




