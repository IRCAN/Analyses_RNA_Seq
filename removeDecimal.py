import sys

mon_fichier = open(sys.argv[1], "r")
contenu = mon_fichier.readlines()

print(str(contenu[0][:-1]))
for lignes in contenu[1:]:
	l=lignes.split()
	a=l[0].split(".")
	print(str(a[0])+"	"+str("	".join(l[1:])))
