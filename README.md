# CART-riboseq-phylogenetics
Github repository for code used to analyze ribo-seq data and build the ribosome rescue factor phylogeny in Tetreault, Callan, Prince, and Feaga, "Collided ribosomes are rescued by the endonuclease Rae1 and trans-translation."

Data Availability

Ribo-seq and RNA-seq data are available on the NCBI Gene Expression Omnibus (GEO) under accession numbers GSE342018 and GSE342019 respectively. csv files including global translatome RPFs (reported in TPM) and RPFs ≥70 nt, as well as ribosome occupancy scores for RPFs ≥70 nt can be found on GEO (GSE342018). Genome accession numbers and presence or absence of each rescue factor, gene sequences used to build HMMER profiles, random genomes in each of the 21 bacterial phyla used to build a maximum-likelihood tree, and all R code used to analyze Ribo-seq data, detect ribosome rescue factor genes, and build the phylogeny are available in this repository. 

csv file contents:

ribosomerescuefactorpresenceandgenomeassemblies.csv:

Accession numbers of genomes that were searched for the presence or absence of ssrA, smpB, rqcH, mutS2, hrpA, smrB, and rae1. Gene presence is indicated as TRUE.

allquerysequencesforrescuefactorphylogeny.csv:

Gene sequences of ssrA, smpB, rqcH, mutS2, hrpA, smrB, and rae1 that were used to build HMMER profiles.

randomgenomes.csv:

Random genomes in each of the 21 bacterial phyla that were used to build a maximum-likelihood tree.



 
