#!/usr/bin/env python3
# coding: utf8

"""
Description:
Using fasta files (scaffold/chromosme/contig file, protein file), gff file, annotation tsv file and the species name
this script writes a genbank file with EC number and Go annotations.

The annotation tsv file contains association between gene and annotation (EC number, GO term)
to add information to the genbank.

The species name needs to be compatible with the taxonomy of the EBI.

Informations need a good formating:
gene ID should be correctly written (like XXX_001 and no XXX_1 if you got more thant 100 genes).
Currently when there is multiple GO terms/EC the script split them when they are separated by ";" or by "," like GO:0006979;GO:0020037;GO:0004601,
if you use another separator add to the re.split(',|;').
For the gff file ensure that the element start position is at least 1.
If it's 0 gffutils will return an error (source : https://github.com/daler/gffutils/issues/104).

Other informations can be added by adding a dictionary with gene ID as key and the information
as value and adapt the condition used for the others annotations (EC, Go term).

"""

import argparse
import datetime
import gffutils
import numpy as np
import os
import pandas as pa
import re
import shutil
import sys

from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from emapper2gbk.utils import is_valid_file, create_GO_namespaces_alternatives, read_annotation, create_taxonomic_data
from typing import Union


def contig_info(contig_id, contig_seq, species_informations):
    """
    Create contig information from species_informations dictionary and contig id and contig seq.
    """
    record = SeqRecord(contig_seq, id=contig_id, name=contig_id,
                    description=species_informations['description'])

    record.seq.alphabet = IUPAC.ambiguous_dna
    if 'data_file_division' in species_informations:
        record.annotations['data_file_division'] = species_informations['data_file_division']
    record.annotations['date'] = datetime.date.today().strftime('%d-%b-%Y').upper()
    if 'topology' in species_informations:
        record.annotations['topology'] = species_informations['topology']
    record.annotations['accessions'] = contig_id
    if 'organism' in species_informations:
        record.annotations['organism'] = species_informations['organism']
    # Use of literal_eval for taxonomy and keywords to retrieve list.
    if 'taxonomy' in species_informations:
        record.annotations['taxonomy'] = species_informations['taxonomy']
    if 'keywords' in species_informations:
        record.annotations['keywords'] = species_informations['keywords']
    if 'source' in species_informations:
        record.annotations['source'] = species_informations['source']

    new_feature_source = sf.SeqFeature(sf.FeatureLocation(1-1,
                                                        len(contig_seq)),
                                                        type="source")
    new_feature_source.qualifiers['scaffold'] = contig_id
    if 'isolate' in species_informations:
        new_feature_source.qualifiers['isolate'] = species_informations['isolate']
    # db_xref corresponds to the taxon NCBI ID.
    # Important if you want to use Pathway Tools after.
    if 'db_xref' in species_informations:
        new_feature_source.qualifiers['db_xref'] = species_informations['db_xref']
    if 'cell_type' in species_informations:
        new_feature_source.qualifiers['cell_type'] = species_informations['cell_type']
    if 'dev_stage' in species_informations:
        new_feature_source.qualifiers['dev_stage'] = species_informations['dev_stage']
    if 'mol_type' in species_informations:
        new_feature_source.qualifiers['mol_type'] = species_informations['mol_type']

    record.features.append(new_feature_source)

    return record


def strand_change(input_strand):
    """
    The input is strand in str ('-', '+') modify it to be a strand in int (-1, +1) to 
    be compatible with SeqIO strand reading.
    """
    if isinstance(input_strand, str):
        if input_strand == '-':
            new_strand = -1
        elif input_strand == '+':
            new_strand = +1
        if input_strand == '.':
            new_strand = None
        elif input_strand == '?':
            new_strand = 0
    elif isinstance(input_strand, int):
        if input_strand == -1:
            new_strand = input_strand
        elif input_strand == +1:
            new_strand = input_strand

    return new_strand


def gff_to_gbk(genome_fasta:str, prot_fasta:str, annotation_data:Union[str, dict], gff_file:str, species_name:str, gbk_out:str, gobasic:Union[None, str, dict]):
    """
    From a genome fasta (containing each contigs of the genome),
    a protein fasta (containing each protein sequence),
    an annotation table (containing gene name associated with GO terms, InterPro and EC),
    a gff file (containing gene, exon, mRNA, ncRNA, tRNA),
    a contig information table (containing species name, taxon ID, ..)
    create a genbank file.
    """

    print('Creating GFF database (gffutils)')
    # Create the gff database file.
    # gffutils use sqlite3 file-based database to access data inside GFF.
    # ':memory:' ask gffutils to keep database in memory instead of writting in a file.
    gff_database = gffutils.create_db(gff_file, ':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

    print('Formatting fasta and annotation file')
    # Dictionary with scaffold/chromosome id as key and sequence as value.
    contig_seqs = OrderedDict()

    for record in SeqIO.parse(genome_fasta, "fasta"):
        id_contig = record.id
        contig_seqs[id_contig] = record.seq


    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seq = {}

    for record in SeqIO.parse(prot_fasta, "fasta"):
        gene_protein_seq[record.id] = record.seq

    # Create a taxonomy dictionary querying the EBI.
    species_informations = create_taxonomic_data(species_name)

    # Read the eggnog tsv file containing GO terms and EC associated with gene name.
    # if metagenomic mode, annotation is already read and given as a dict
    if not type(annotation_data) is dict:
        annotation_data = dict(read_annotation(annotation_data))

    # Query Gene Ontology to extract namespaces and alternative IDs.
    # go_namespaces: Dictionary GO id as term and GO namespace as value.
    # go_alternatives: Dictionary GO id as term and GO alternatives id as value.
    if not type(gobasic[0]) is dict and not type(gobasic[1]) is dict:
        go_namespaces, go_alternatives = create_GO_namespaces_alternatives(gobasic)
    else:
        go_namespaces, go_alternatives = gobasic

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    seq_objects = []

    print('Assembling Genbank informations')

    # Iterate through each contig.
    # Then iterate through gene and throug RNA linked with the gene.
    # Then look if protein informations are available.
    for gene in gff_database.features_of_type('gene'):
        id_gene = gene.id
        if id_gene.isnumeric():
            id_gene = f"gene_{id_gene}"
        elif "|" in id_gene:
            id_gene = id_gene.split("|")[0]
        else:
            id_gene = id_gene

        record = contig_info(id_gene, contig_seqs[id_gene], species_informations)

        start_position = gene.start -1
        end_position = gene.end
        strand = strand_change(gene.strand)
        new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                            type="gene")
        new_feature_gene.qualifiers['locus_tag'] = id_gene
        # Add gene information to contig record.
        record.features.append(new_feature_gene)

        # Iterate through gene childs to find CDS object.
        # For each CDS in the GFF add a CDS in the genbank.
        for cds_object in gff_database.children(gene, featuretype="CDS", order_by='start'):
            cds_id = cds_object.id
            start_position = cds_object.start -1
            end_position = cds_object.end
            new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                end_position,
                                                                strand),
                                                            type="CDS")

            new_feature_cds.qualifiers['locus_tag'] = id_gene

            # Add GO annotation according to the namespace.
            if id_gene in annotation_data.keys():
                # Add gene name.
                if 'Preferred_name' in annotation_data[id_gene]:
                    new_feature_cds.qualifiers['gene'] = annotation_data[id_gene]['Preferred_name']

                if 'GOs' in annotation_data[id_gene] :
                    gene_gos = annotation_data[id_gene]['GOs'].split(',')
                    if gene_gos != [""]:
                        go_components = []
                        go_functions = []
                        go_process = []

                        for go in gene_gos:
                            # Check if GO term is not a deprecated one.
                            # If yes take the corresponding one in alternative GO.
                            if go not in go_namespaces:
                                go_test = go_alternatives[go]
                            else:
                                go_test = go
                            if go_namespaces[go_test] == 'cellular_component':
                                    go_components.append(go)
                            if go_namespaces[go_test] == 'molecular_function':
                                go_functions.append(go)
                            if go_namespaces[go_test] == 'biological_process':
                                go_process.append(go)
                        new_feature_cds.qualifiers['go_component'] = go_components
                        new_feature_cds.qualifiers['go_function'] = go_functions
                        new_feature_cds.qualifiers['go_process'] = go_process


                # Add EC annotation.
                if 'EC' in annotation_data[id_gene]:
                    gene_ecs = annotation_data[id_gene]['EC'].split(',')
                    if gene_ecs != [""]:
                        new_feature_cds.qualifiers['EC_number'] = gene_ecs

            if id_gene in gene_protein_seq:
                new_feature_cds.qualifiers['translation'] = gene_protein_seq[id_gene]

            # Add CDS information to contig record
            record.features.append(new_feature_cds)

        seq_objects.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(seq_objects, gbk_out, 'genbank')

def main(genome_fasta, prot_fasta, annot_table, gff_file, species_name, gbk_out, gobasic=None):
    # check validity of inputs
    for elem in [genome_fasta, prot_fasta]:
        print(elem)
        if not is_valid_file(elem):
            print(f"{elem} is not a valid path file.")
            sys.exit(1)

    gff_to_gbk(genome_fasta, prot_fasta, annot_table, gff_file, species_name, gbk_out, gobasic)
