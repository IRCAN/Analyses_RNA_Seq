# PipelineNGS

I- Liste des fichiers

Analyse_kmer.R
config.txt
differential_expression.R
legend_image.py
PipelineNGS.py
README.md
recup_info.py
recup_nom_gene.R
removeDecimal.py
sleuth.R

II- Description

Ce programme permet de réaliser la plupart des étapes nécessaires à l'analyses de données de RNA-seq.

Mapping avec le logiciel STAR
Comptage des reads par gènes avec Featurecounts
Analyse de l'expression différentielle avec le package bioconductor DESeq2
Enrichissement, détection de pathways


Toutes les étapes peuvent être réalisées séparement. Néanmoins, il est recommandé de réaliser les étapes dans l'ordre chronologique, de manière à ce que les fichiers de sortie d'un étape n correspondent aux fichiers d'entrées pour l'étape n+1.

III-Ressources nécessaires

Avant de lancer le pipeline, les fichiers ci-dessous sont nécessaires pour son bon fonctionnement.

-Avoir réalise run index STAR du génome ou transcriptome sur lequel on veut réaliser le mapping. Si pipeline lancé du serveur de l'ircan, les index pour l'humain (Hg19- Hg38),pour la souris (Mg38) et pour le poissson-zèbre (zv10) sont déjà disponibles.
Si extérieur IRCAN, ou pour tout nouveau génome, il faudra modifier les chemins dans le fichier config.txt

-le fichier .gtf correspondant, à récupérer sur GENCODE (déjà disponibles pour IRCAN); chemin également à modifier dans le fichier config.txt

-Pour l'expression différentielle, un fichier config.txt doit être fourni avant le lancement du pipeline.
Le fichier doit être comme ci-dessous

----------------------------------------------------------------------------------------------------------

Tell if reads are paired-end (if Paired-end=FALSE -> reads = single-end):
#Do not put space after "="
>Paired-end =TRUE

Choose your genome:
#Do not put space after "="
>Genome =mm10

Available genomes:
Human:	hg37,hg38	Mouse: mm10,mm9 poissonzebre: zv10
		
############################################################################
#############                                 ##############################
#############      Differential Expression    ##############################
#############                                 ##############################
############################################################################

Expression différentielle =True
--------Description échantillons---------------(put tab between sample and condition)
sample	condition
NameSample6	C1
NameSample1	C1
NameSample3	C2
NameSample4	C2
NameSample2	C3
NameSample5	C3
-----------conditions à tester-------------------(put tab between the two conditions)
condition1	condition2
C1	C2
C2	C3
----------------------------------------------------------------------------------------------------------------------------------------------
III-


Étapes (ordre chronologique)

a- Création des répertoires pour stocker les fichiers de sorties de chaque étape.
Les noms des répertoires propres à chaque échantillons sont basés sur les mêmes noms que ceux de leur fichier fastq respectif.


b- Analyses des kmers (supprimé?)

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


Avec le logiciel deseq2, une analyse de l'expression différentielle est effectuée entre 2 conditions (selon fichier config.txt)
la liste des gènes différenitellement exprimées est ensuite enrichie avec les identifiants Kegg, les GO terms, ...
Ensuite les Kegg pathways avec le plus de gènes D.E. impliqués sont créés.


IV- Lancement rapide

Pour lancer le programme avec toutes les étapes, lancer la commande ci-dessous dans le terminal

python PipelineNGS.py -i DATA -o OUTPUT -s all -c config.txt




**Requierement**

-STAR (Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29:15–21. )

-featurecounts (Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30, 2014)

-Salmon (Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.)

-DESeq2  (Love MI, Huber W and Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, pp. 550. doi: 10.1186/s13059-014-0550-8. )

-Gage (Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161)
-Pathwiew  (Weijun Luo and Cory Brouwer. Pathview: an R/Bioconductor package for pathway-based data integration and visualization. Bioinformatics, 29(14):1830-1831, 2013. doi: 10.1093/bioinformatics/btt285 )

# PipelineNGS
