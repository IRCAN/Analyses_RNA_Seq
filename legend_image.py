#python LegendImage.py mmu04974liste_gene.txt mmu04974.Upregulate_1PV_.png keggres.txt 


from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
import sys



def recup_legend_image(liste_gene,img_input, keggres):
	numero_path=img_input.split(".")[0]


	mon_fichier = open(liste_gene, "r")
	contenu = mon_fichier.readlines()

	listeGeneName=[]
	n=contenu[0][:-1].split("\t")
	try:
		rang_nom_gen=n.index("\"gene_name\"")
		rang_path=n.index("\"PATH\"")

	except:
		rang_nom_gen=n.index("gene_name")
		rang_path=n.index("PATH")

	

	for lignes in contenu[1:]:
		lignessplit=lignes.split("\t")
		try:
			if numero_path[3:] in lignessplit[rang_path]:	
				nameGene=lignessplit[rang_nom_gen]
				listeGeneName.append(nameGene)
		except:
			pass
	
	print(listeGeneName)
	a=str(', '.join(listeGeneName))
	i=0
	b=""
	for caract in a:
		i=i+1
		if caract==" ":
			if i>150:
				i=0
				b=b+"\n"
		else:
			b=b+caract
	

	mon_fichier2 = open(keggres, "r")
	contenu2 = mon_fichier2.readlines()
	qvalue=""
	for lignes in contenu2[1:]:
		lignessplit=lignes.split("\t")
		refPathway=str(lignessplit[0][1:9])

		if refPathway in str(numero_path):

			qvalue=format(float(lignessplit[4]), '.3f')
		

	ma_chaine = "Genes involved in pathway: "+str(len(listeGeneName))+" \n \n"+str(b)+" \n \n"+"q-value="+qvalue

	image_in = Image.open(img_input)
	w, h = image_in.size

	my_image = Image.new("RGBA", (w, h+200), "white")

	my_image.paste(image_in, (0,0))

	carac = ImageDraw.Draw(my_image)
	#font = ImageFont.truetype("arial.ttf", 20)
	carac.text( (10, h+20), unicode(ma_chaine,'UTF-8'),  fill="black")
	im_output="Qvalue_"+qvalue+"_"+str(numero_path)+".png"
	my_image.save(im_output,"png")



recup_legend_image(sys.argv[1],sys.argv[2],sys.argv[3])
