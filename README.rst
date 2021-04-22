.. image:: https://img.shields.io/pypi/v/emapper2gbk.svg
	:target: https://pypi.org/project/emapper2gbk

.. image:: https://img.shields.io/github/license/AuReMe/emapper_to_gbk.svg
	:target: https://github.com/AuReMe/emapper_to_gbk/blob/master/LICENSE

.. image:: https://github.com/AuReMe/emapper_to_gbk/workflows/Python%20package/badge.svg
    :target: https://github.com/AuReMe/emapper_to_gbk/actions

.. image:: https://img.shields.io/badge/doi-10.7554/eLife.61968-blueviolet.svg
	:target: https://doi.org/10.7554/eLife.61968

emapper2gbk: creation of genbank files from Eggnog-mapper annotation outputs
============================================================================

Starting from fasta and `Eggnog-mapper <http://eggnog-mapper.embl.de/>`__ annotation files, build a gbk file that is suitable for metabolic network reconstruction with `Pathway Tools <http://bioinformatics.ai.sri.com/ptools/>`__. Adds the GO terms and EC numbers annotations in the genbank file.
There are two main modes:

* **genes mode**: suitable when a list of isolated genes/proteins have been annotated with Eggnog-mapper, typically the gene catalogue of a metagenome.

* **genomes mode**: usually when focusing on a single organism, with a ``.gff`` file. The creation of genbanks can be performed in parallel by providing directories (with matching names for genomes, proteomes and annotation files) as inputs.

**If you use emapper2gbk, please cite**

Belcour* A, Frioux* C, Aite M, Bretaudeau A, Hildebrand F, Siegel A. Metage2Metabo, microbiota-scale metabolic complementarity for the identification of key species. eLife 2020;9:e61968 `https://doi.org/10.7554/eLife.61968 <https://doi.org/10.7554/eLife.61968>`_ .

Main inputs
-----------

emapper2gbk genes
~~~~~~~~~~~~~~~~~

For each annotated list of genes, inputs are:

* a nucleotide fasta file containing the CDS sequence of each genes or a folder containing multiple nucleotide fasta files.
* the translated sequences in amino-acids in fasta or a folder containing the corresponding protein sequences to the nucleotide sequences (must be the same name).
* the annotation file obtained after Eggnog-mapper annotation (usually ``xxx.emapper.annotation``) or a folder with multiple annotation files (must be the same name as nucleotide fasta file and ends with '.tsv' extension).

In addition, as optional files:

* the name of the considered organism (can be "bacteria" or "metagenome") or a file with organisms names (matching the genomes names).
* the merge option to merge genes into fake contigs.
* the number of available cores for multiprocessing (when working on multiple genomes).
* a go-basic file of GO ontology (if not given, emapper2gbk will download a copy and use it).

Example:

Input with files:

.. code-block:: text

    nucleotide_sequences.fna
    protein_sequence.faa
    annotation.emapper.annotation

Input with folders:

.. code-block:: text

    nucleotide_sequences
    ├── gene_list_1.fna
    ├── gene_list_2.fna
    protein_sequence
    ├── gene_list_1.faa
    ├── gene_list_2.faa
    annotation
    ├── gene_list_1.tsv
    ├── gene_list_2.tsv

.. image:: pictures/emapper2gbk_genes.svg

To work the ID of the genes in the nucleic fasta file (``-fn``) must be the same than the ID of the proteins in the protein fasta file (``-fp``) and in the annotation file (``-a``).

emapper2gbk genomes
~~~~~~~~~~~~~~~~~~~

For each genomes, inputs are:

* a nucleotide fasta file containing the sequence of each contigs/chromosomes for the genome or a folder containing multiple nucleotide fasta files.
* the proteome corresponding to the genome or a folder containing the corresponding protein sequences to the nucleotide sequences (must be the same name).
* the GFF file corresponding to the genome or a folder containing multiple GFF files (must be the same name as the nucleotide folder).
* the annotation file obtained after Eggnog-mapper annotation (usually ``xxx.emapper.annotation``) or a folder with multiple annotation files (must be the same name as nucleotide fasta file and ends with '.tsv' extension)

In addition, as optional files:

* the name of the considered organism (can be "bacteria") or a file with organisms names (matching the genomes names).
* the number of available cores for multiprocessing (when working on multiple genomes).
* a go-basic file of GO ontology (if not given, emapper2gbk will download a copy and use it).

Example:

Input with files:

.. code-block:: text

    nucleotide_sequences.fna
    protein_sequence.faa
    annotation.emapper.annotation
    genome.gff

Input with folders:

.. code-block:: text

    nucleotide_sequences
    ├── genome_1.fna
    ├── genome_2.fna
    protein_sequence
    ├── genome_1.faa
    ├── genome_2.faa
    annotation
    ├── genome_1.tsv
    ├── genome_2.tsv
    gff
    ├── genome_1.gff
    ├── genome_2.gff

.. image:: pictures/emapper2gbk_genomes.svg

The ID in the chromosome/contigs/scaffolds fasta file (``-fn``) must correspond to region in the gff file (``-g``).
Then the genes in the region will be found and the child CDS associated to the genes wil be extracted.
The CDS ID must be the same than the ID in the protein fasta file (``-fp``) and the ID in the eggnog-mapper annotation file (``-a``).

Dependencies and installation
-----------------------------

Dependencies
~~~~~~~~~~~~

All are described in ``requirements.txt`` and can be installed with ``pip install -r requirements.txt``.

* biopython
* gffutils
* pandas
* pronto
* requests

Install
~~~~~~~

* From this cloned repository

.. code-block:: sh

    pip install -r requirements.txt
    pip install .

* From Pypi

.. code-block:: sh

    pip install emapper2gbk

Usage
-----

Convert GFF, fastas, annotation table and species name into Genbank.

.. code-block:: sh

    usage: emapper2gbk [-h] [-v] {genes,genomes} ...

    Starting from fasta and Eggnog-mapper annotation files, build a gbk file that is suitable for metabolic network reconstruction with Pathway Tools. Adds the GO terms and EC numbers annotations in the genbank file.

    Two modes:
    - genomes (one genome/proteome/gff/annot file --> one gbk).
    - genes with the annotation of the full gene catalogue and fasta files (nucleic and protein) corresponding to list of genes.

    Examples:

    * Genomic - single mode

    emapper2gbk genomes -fn genome.fna -fp proteome.faa -gff genome.gff -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo]

    * Genomic - multiple mode, "bacteria" as default name

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -n metagenome -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

    * Genomic - multiple mode, tsv file for organism names

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -nf matching_genome_orgnames.tsv -o gbk_dir/ -a eggnog_annotation_dir/ [-go go-basic.obo]

    * Metagenomic

    emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv --one-annot-file [-go go-basic.obo]

    You can give the GO ontology as an input to the program, it will be otherwise downloaded during the run. You can download it here: http://purl.obolibrary.org/obo/go/go-basic.obo .
    The program requests the NCBI database to retrieve taxonomic information of the organism. However, if the organism is "bacteria" or "metagenome", the taxonomic information will not have to be retrieved online.
    Hence, if you need to run the program from a cluster with no internet access, it is possible for a "bacteria" or "metagenome" organism, and by providing the GO-basic.obo file.
    For specific help on each subcommand use: emapper2gbk {cmd} --help

    optional arguments:
    -h, --help       show this help message and exit
    -v, --version    show program's version number and exit

    subcommands:
    valid subcommands:

    {genes,genomes}
        genes          genomic mode : 1-n annot, 1-n faa, 1-n fna (gene sequences) --> 1 gbk
        genomes        genoems mode: 1 contig/chromosome fasta, 1 protein fasta, 1 GFF, 1 annot --> 1 gbk


* Genomes mode

  * Usage

    .. code-block:: sh

        usage: emapper2gbk genomes [-h] -fn FASTANUCLEIC -fp FASTAPROT -o OUPUT_DIR [-g GFF] [-nf NAMEFILE] [-n NAME] -a
                                ANNOTATION [-c CPU] [-go GOBASIC] [-q]

        Use the annotation of a complete gene catalogue and build gbk files for each set of genes (fna) and proteins (faa)
        from input directories

        optional arguments:
        -h, --help            show this help message and exit
        -fn FASTANUCLEIC, --fastanucleic FASTANUCLEIC
                                fna file or directory
        -fp FASTAPROT, --fastaprot FASTAPROT
                                faa file or directory
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory/file path
        -g GFF, --gff GFF     gff file or directory
        -nf NAMEFILE, --namefile NAMEFILE
                                organism/genome name (col 2) associated to genome file basenames (col 1). Default =
                                'metagenome' for metagenomic and 'cellular organisms' for genomic
        -n NAME, --name NAME  organism/genome name in quotes
        -a ANNOTATION, --annotation ANNOTATION
                                eggnog annotation file or directory
        -c CPU, --cpu CPU     cpu number for metagenomic mode or genome mode using input directories
        -go GOBASIC, --gobasic GOBASIC
                                go ontology, GOBASIC is either the name of an existing file containing the GO Ontology or the
                                name of the file that will be created by emapper2gbk containing the GO Ontology
        -q, --quiet           quiet mode, only warning, errors logged into console

  * Examples

    * Genomic - single mode

    .. code:: sh

      emapper2gbk genomes -fn genome.fna -fp proteome.faa -gff genome.gff -n "Escherichia coli" -o coli.gbk -a eggnog_annotation.tsv [-go go-basic.obo]

    * Genomic - multiple mode, "bacteria" as default name

* genes mode

  * Usage

    .. code-block:: sh

        usage: emapper2gbk genes [-h] -fn FASTANUCLEIC -fp FASTAPROT -o OUPUT_DIR [--one-annot-file] -a ANNOTATION [-c CPU]
                                [-n NAME] [-nf NAMEFILE] [-go GOBASIC] [--merge MERGE] [-q]

        Build a gbk file for each genome/set of genes with an annotation file for each

        optional arguments:
        -h, --help            show this help message and exit
        -fn FASTANUCLEIC, --fastanucleic FASTANUCLEIC
                                fna file or directory
        -fp FASTAPROT, --fastaprot FASTAPROT
                                faa file or directory
        -o OUPUT_DIR, --out OUPUT_DIR
                                output directory/file path
        --one-annot-file      Option to use when there is only one annotation file for multiples genes fastas.
        -a ANNOTATION, --annotation ANNOTATION
                                eggnog annotation file or directory
        -c CPU, --cpu CPU     cpu number for metagenomic mode or genome mode using input directories
        -n NAME, --name NAME  organism/genome name in quotes
        -nf NAMEFILE, --namefile NAMEFILE
                                organism/genome name (col 2) associated to genome file basenames (col 1). Default =
                                'metagenome' for metagenomic and 'cellular organisms' for genomic
        -go GOBASIC, --gobasic GOBASIC
                                go ontology, GOBASIC is either the name of an existing file containing the GO Ontology or the
                                name of the file that will be created by emapper2gbk containing the GO Ontology
        --merge MERGE         Number of gene sequences to merge into fake contig from a same file in the genbank file.
        -q, --quiet           quiet mode, only warning, errors logged into console

  * Example

    .. code:: sh

      emapper2gbk genes -fn genome_dir/ -fp proteome_dir/ -o gbk_dir/ -a gene_cat_ggnog_annotation.tsv [-go go-basic.obo]

