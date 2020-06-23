#!/usr/bin/env python3
# coding: utf8
import sys
import argparse
import os
import re
import shutil
import logging
from typing import Union
from collections import OrderedDict
from Bio import SeqFeature as sf
from Bio import SeqIO

from emapper2gbk.utils import is_valid_file, create_GO_namespaces_alternatives, read_annotation, create_taxonomic_data, get_basename, record_info, create_cds_feature

logger = logging.getLogger(__name__)

"""
Description:
Using fasta files (scaffold/chromosme/contig file, protein file), annotation tsv file from eggnog and the species name
this script writes a genbank file with EC number and Go annotations.
The species name needs to be compatible with the taxonomy of the EBI.
Informations need a good formating:
gene ID should be correctly written (like XXX_001 and no XXX_1 if you got more thant 100 genes).
Currently when there is multiple GO terms/EC the script split them when they are separated by ";" or by "," like GO:0006979;GO:0020037;GO:0004601,
if you use another separator add to the re.split(',|;').
Other informations can be added by adding a dictionary with gene ID as key and the information
as value and adapt the condition used for the others annotations (EC, Go term).
"""


def faa_to_gbk(genome_fasta:str, prot_fasta:str, annotation_data:Union[str, dict], species_name:str, gbk_out:str, gobasic:Union[None, str, dict]):
    """
    From a genome fasta (containing each genes of the genome),
    a protein fasta (containing each protein sequence),
    an eggnog annotation table,
    a species name,
    a name for the genbank output,
    gobasic: name of go-basic.obo or a tuple of dictionaries (go_namespaces, go_alternatives),
    create a genbank file.
    """
    genome_id = get_basename(genome_fasta)

    logger.info('Formatting fasta and annotation file for ' + genome_id)

    # Dictionary with gene id as key and nucleic sequence as value.
    gene_nucleic_seqs = OrderedDict()

    for record in SeqIO.parse(genome_fasta, "fasta"):
        gene_nucleic_seqs[record.id] = record.seq

    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seqs = SeqIO.to_dict(SeqIO.parse(prot_fasta, 'fasta'))

    for record in SeqIO.parse(prot_fasta, "fasta"):
        gene_protein_seqs[record.id] = record.seq

    # Create a taxonomy dictionary querying the EBI.
    species_informations = create_taxonomic_data(species_name)

    # Read the eggnog tsv file containing GO terms and EC associated with gene name.
    # if metagenomic mode, annotation is already read and given as a dict
    if not type(annotation_data) is dict:
        annotation_data = dict(read_annotation(annotation_data))

    # Query Gene Ontology to extract namespaces and alternative IDs.
    # go_namespaces: Dictionary GO id as term and GO namespace as value.
    # go_alternatives: Dictionary GO id as term and GO alternatives id as value.
    if gobasic:
        if not type(gobasic[0]) is dict and not type(gobasic[1]) is dict:
            go_namespaces, go_alternatives = create_GO_namespaces_alternatives(gobasic)
        else:
            go_namespaces, go_alternatives = gobasic
    else:
        go_namespaces, go_alternatives = create_GO_namespaces_alternatives()

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    seq_objects = []

    logger.info('Assembling Genbank informations for ' + genome_id)

    # Iterate through each contig/gene.
    for gene_nucleic_id in sorted(gene_nucleic_seqs):
        # Create a SeqRecord object using gene information.
        record = record_info(gene_nucleic_id, gene_nucleic_seqs[gene_nucleic_id], species_informations)
        # if id is numeric, change it
        if gene_nucleic_id.isnumeric():
            id_gene = f"gene_{gene_nucleic_id}"
        elif "|" in gene_nucleic_id:
            id_gene = gene_nucleic_id.split("|")[0]
        else:
            id_gene = gene_nucleic_id
        start_position = 1
        end_position = len(gene_nucleic_seqs[gene_nucleic_id])
        strand = 0
        new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                            type="gene")
        new_feature_gene.qualifiers['locus_tag'] = id_gene  + "_0001"
        # print(new_feature_gene.qualifiers['locus_tag'] )
        # Add gene information to contig record.
        record.features.append(new_feature_gene)

        new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                            type="CDS")

        # Add gene ID in locus_tag.
        new_feature_cds.qualifiers['locus_tag'] = id_gene + "_0001"

        # Add functional annotations.
        if gene_nucleic_id in annotation_data.keys():
            # Add gene name.
            if 'Preferred_name' in annotation_data[gene_nucleic_id]:
                new_feature_cds.qualifiers['gene'] = annotation_data[gene_nucleic_id]['Preferred_name']

            # Add GO annotation according to the namespace.
            if 'GOs' in annotation_data[gene_nucleic_id]:
                gene_gos = annotation_data[gene_nucleic_id]['GOs'].split(',')
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
            if 'EC' in annotation_data[gene_nucleic_id]:
                gene_ecs = annotation_data[gene_nucleic_id]['EC'].split(',')
                if gene_ecs != [""]:
                    new_feature_cds.qualifiers['EC_number'] = gene_ecs

        if gene_nucleic_id in gene_protein_seqs.keys():
            # Add protein sequence.
            new_feature_cds.qualifiers['translation'] = gene_protein_seqs[gene_nucleic_id]

        # Add CDS information to contig record
        record.features.append(new_feature_cds)

        seq_objects.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(seq_objects, gbk_out, 'genbank')


def main(genome_fasta, prot_fasta, annot_table, species_name, gbk_out, gobasic=None):
    # check validity of inputs
    for elem in [genome_fasta, prot_fasta]:
        if not is_valid_file(elem):
            logger.critical(f"{elem} is not a valid path file.")
            sys.exit(1)

    faa_to_gbk(genome_fasta, prot_fasta, annot_table, species_name, gbk_out, gobasic)
