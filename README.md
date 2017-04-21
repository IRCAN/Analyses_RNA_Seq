# PipelineNGS

I- Liste des fichiers

Analyse_kmer.R
differential_expression.R	
PipelineNGS.py
README.md
Recup_info.py


II- Description

Ce programme permet de réaliser la plupart des étapes nécessaires à l'analyses de données de RNA-seq.

Mapping avec le logiciel STAR
Comptage des reads par gènes avec Featurecounts
Analyse de l'expression différentielle avec le package bioconductor DESeq2
Enrichissement, détection de pathways


Toutes les étapes peuvent être réalisées séparement. Néanmoins, il est recommandé de réaliser les étapes dans l'ordre chronologique, de manière à ce que les fichiers de sortie d'un étape n correspondent aux fichiers d'entrées pour l'étape n+1.

III-Ressources nécessaires

Avant de lancer le pipeline, les fichiers ci-dessous sont nécessaires pour son bon fonctionnement.

-Avoir réalise run index STAR du génome ou transcriptome sur lequel on veut réaliser le mapping. Si pipeline lancé du serveur de l'ircan, les index pour l'humain (Hg19- Hg38) et pour la souris (Mg38) sont déjà disponibles.
Si extérieur IRCAN, ou pour tout nouveau génome, il faudra modifier les chemins dans le script PipelineNGS.py

-le fichier .gtf correspondant, à récupérer sur GENCODE (déjà disponibles pour IRCAN)

-Pour l'expression différentielle, un fichier serie.txt doit être fourni avant le lancement du pipeline.
Le fichier doit être comme ci-dessous, avec le sampleName qui doivent correspondrent aux noms des échantillons donnés pour les fichiers fastq.

sampleName	condition	replicat
NameSample6	C1	R1
NameSample1	C1	R2
NameSample3	C2	R1
NameSample4	C2	R1
NameSample2	C3	R1
NameSample5	C3	R2

III-


Étapes (ordre chronologique)

a- Création des répertoires pour stocker les fichiers de sorties de chaque étape.
Les noms des répertoires propres à chaque échantillons sont basés sur les mêmes noms que ceux de leur fichier fastq respectif.


b- Analyses des kmers

Un premier arbre de classification des données est réalisé, à partir des outils de la librairie SeqTools (R).
Cet arbre se fait par l'analyse des kmers des fichiers fastq, et ne nécessite donc pas un mapping des échantillons au préalable.
Indépendant, les fichiers de sorties ne seront pas réutilisés pour une autre étape.

INPUT:
répertoire contenant les fichiers fastq
OUTPUT:
Rplots.pdf : arbre des échantillons


c-Alignement

Alignement des échantillons par le logiciel STAR 
les fichiers .Bam sont placés dans les répertoires de chaque échantillon

INPUT:
Répertoire contenant les fichiers fastq des échantillons
OUTPUT:
Un répertoire pour chaque échantillon avec tous les fichiers de sortie fournis par STAR

d- featurecounts
INPUT: fichiers .bam pour chaque échantillons

OUTPUT:
un tableu de comptage de reads pour chaque gene_id pour chaque échantillon, présent dans chaque répertoire
Un tableau regroupant tous les différents comptages: allcounts.txt 


e-Expression différentielle, Annotations




IV- Lancement rapide

Pour lancer le programme avec toutes les étapes, lancer la commande ci-dessous dans le terminal

python PipelineNGS.py /chemin/répertoire/contenant/fastq/ /chemin/nom/répertoire/fichiers/sorties -id [Mouse ou Human] -seq [single ou paired] --All ALL






# PipelineNGS
