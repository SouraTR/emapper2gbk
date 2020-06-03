# Creation of genbank files from Eggnog-mapper annotation

Starting from fasta and [Eggnog-mapper](http://eggnog-mapper.embl.de/) annotation files, build a gbk file that is suitable for metabolic network reconstruction with Pathway Tools. Adds the GO terms and EC numbers annotations in the genbank file.
The program can be run with or without a `.gff` file.
There are two main modes:

* **genomic mode**: usually when focusing on a single organism, with a `.gff` file (but can be without). The creation of genbanks can be performed in parallel by providing directories (with matching names for genomes, proteomes and annotation files) as inputs

* **metagenomic mode**: suitable when a list of isolated genes/proteins have been annotated with Eggnog-mapper, typically the gene catalogue of a metagenome. In that case, there is one annotation file and multiple genomes/proteomes that will contain a subset of the genes present in the annotated gene catalogue. The purpose is to create a genbank file with the genes and annotations that matches each genome (Metagenomic species core genome, Metagenomic-assembled genome...) content. There will not be any `.gff` associated to the genomes, and the creation of genbanks can be performed in parallel by providing directories (with matching names for genome and proteome files) as inputs. In that case, you can give as organism name "metagenome" or "bacteria". 

## Main inputs

For each annotated genome, inputs consist of (but are not limited to):
* a nucleotide fasta file
* the translated sequences in fasta
* the annotation file obtained after Eggnog-mapper annotation (usually `xxx.emapper.annotation`)
* the name of the considered organism (can be "bacteria" or "metagenome") or a file with organisms names (matching the genomes names).

In addition, the gff file, the number of available cores for multiprocessing, the go-basic file of GO ontology

## Dependencies and installation


## Usage 

Convert GFF, fastas, annotation table and species name into Genbank.
usage: emapper2gbk [-h] [-v] {genomic,metagenomic} ...



Two modes: genomic (one genome/proteome/gff/annot file --> one gbk) or metagenomic with the annotation of the full gene catalogue and fasta files (proteome/genomes) corresponding to list of genes.

Genomic mode can be used with or without gff files making it suitable to build a gbk from a list of genes and their annotation.

Examples:

* Genomic - single mode

emapper2gbk genomic -fg genome.fna -fp proteome.faa [-gff genome.gff] -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo]

* Genomic - multiple mode, "bacteria" as default name

emapper2gbk genomic -fg genome_dir/ -fp proteome_dir/ [-gff gff_dir/] -n metagenome -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

* Genomic - multiple mode, tsv file for organism names

emapper2gbk genomic -fg genome_dir/ -fp proteome_dir/ [-gff gff_dir/] -nf matching_genome_orgnames.tsv -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

* Metagenomic

emapper2gbk metagenomic -fg genome_dir/ -fp proteome_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv [-go go-basic.obo]

You can give the GO ontology as an input to the program, it will be otherwise downloaded during the run. You can download it here: http://purl.obolibrary.org/obo/go/go-basic.obo .
The program requests the NCBI database to retrieve taxonomic information of the organism. However, if the organism is "bacteria" or "metagenome", the taxonomic information will not have to be retrieved online.
Hence, if you need to run the program from a cluster with no internet access, it is possible for a "bacteria" or "metagenome" organism, and by providing the GO-basic.obo file.
 For specific help on each subcommand use: emapper2gbk {cmd} --help

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  valid subcommands:

  {genomic,metagenomic}
    genomic             genomic mode : 1 annot, 1 faa, 1 fna, [1 gff] --> 1 gbk
    metagenomic         metagenomic mode : 1 annot, n faa, n fna --> n gb