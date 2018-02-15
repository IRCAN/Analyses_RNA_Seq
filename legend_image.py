#python LegendImage.py mmu04974liste_gene.txt mmu04974.Upregulate_1PV_.png keggres.txt 


from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import sys



import sys
def recup_list_gene(file1, file2):
	titre=file1
	mon_fichier = open(file2, "r")
	contenu = mon_fichier.readlines()

	path2=titre.split("/")[-1]
	path=path2.split(".")[0][3:]

	liste=[]
	b=0
	for lignes in contenu:

		l=lignes.split("\t")
		if path in l[13]:
			liste.append(l[1])
			b=b+1

	return(liste, path)
	





def recup_legend_image(liste_gene,numero_path,img_input, keggres, output):
	
	a=str(', '.join(liste_gene))
	aa= a.replace('"','') 
	mot=aa.split(",")
	mot2=sorted(mot)
	i=0
	b=""
	for m in mot2:
		i=i+len(m)
		if i>150:
			i=0
			b=b+m+"\n"
		else:
			b=b+m
		#if caract==" ":
		#	if i>150:
		#		i=0
		#		b=b+"\n"
		#elif caract!='"':
		#	i=i+1
		#	b=b+caract
	

	mon_fichier2 = open(keggres, "r")
	contenu2 = mon_fichier2.readlines()
	qvalue=""
	for lignes in contenu2[1:]:
		lignessplit=lignes.split("\t")
		refPathway=str(lignessplit[0][1:9])
		if str(numero_path) in refPathway:

			qvalue=format(float(lignessplit[4]), '.3f')
		

	ma_chaine = "Genes involved in pathway: "+str(len(liste_gene))+" \n \n"+str(b)+" \n \n"+"q-value="+qvalue

	image_in = Image.open(img_input)
	w, h = image_in.size

	my_image = Image.new("RGBA", (w, h+200), "white")

	my_image.paste(image_in, (0,0))

	carac = ImageDraw.Draw(my_image)
	#font = ImageFont.truetype("arial.ttf", 20)
	carac.text( (10, h+20), unicode(ma_chaine,'UTF-8'),  fill="black")
	im_output=output+"/Qvalue_"+qvalue+"_"+str(numero_path)+".png"
	my_image.save(im_output,"png")


liste_gene,numero_path=recup_list_gene(sys.argv[1],sys.argv[2])
recup_legend_image(liste_gene,numero_path,sys.argv[1],sys.argv[3], sys.argv[4])
