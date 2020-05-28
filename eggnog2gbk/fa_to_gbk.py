#!/usr/bin/env python3
# coding: utf8
import sys
import argparse
import datetime
import os
import re
import shutil
import logging
from typing import Union
from collections import OrderedDict
from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from eggnog2gbk.utils import is_valid_file, create_GO_dataframes, read_annotation, create_taxonomic_data

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

def contig_info(contig_id, contig_seq, species_informations):
    """
    Create contig information from species_informations dictionary and contig id and contig seq.
    """
    if contig_id.isnumeric():
        newname = f"_{contig_id}"
    elif "|" in contig_id:
        newname = contig_id.split("|")[0]
    else:
        newname = contig_id
    record = SeqRecord(contig_seq, id=contig_id, name=newname,
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

def faa_to_gbk(genome_fasta:str, prot_fasta:str, annotation_data:Union[str, dict], species_name:str, gbk_out:str, gobasic:str=None):
    """
    From a genome fasta (containing each contigs of the genome),
    a protein fasta (containing each protein sequence),
    an annotation table (containing gene name associated with GO terms, InterPro and EC),
    a gff file (containing gene, exon, mRNA, ncRNA, tRNA),
    a contig information table (containing species name, taxon ID, ..)
    create a genbank file.
    """
    logger.info('Formatting fasta and annotation file')
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

    # Read the ggnog tsv file containing GO terms and EC associated with gene name.
    # if metagenomic mode, annotation is already read and given as a dict
    if not type(annotation_data) is dict:
        annotation_data = read_annotation(annotation_data)

    # Query Gene Ontology to extract namespaces and alternative IDs.
    df_go_namespace, df_go_alternative = create_GO_dataframes(gobasic)
    # Dictionary GO id as term and GO namespace as value.
    df_go_namespace.set_index('GO', inplace=True)
    go_namespaces = df_go_namespace['namespace'].to_dict()

    # Dictionary GO id as term and GO alternatives id as value.
    df_go_alternative.set_index('GO', inplace=True)
    go_alternatives = df_go_alternative['alternative_GO'].to_dict()

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    seq_objects = []

    logger.info('Assembling Genbank informations')

    # Iterate through each contig/gene.
    for contig_id in sorted(contig_seqs):
        # Data for each contig.
        record = contig_info(contig_id, contig_seqs[contig_id], species_informations)
        # if id is numeric, change it
        if contig_id.isnumeric():
            id_gene = f"gene_{contig_id}"
        elif "|" in contig_id:
            id_gene = contig_id.split("|")[0]
        else:
            id_gene = contig_id
        start_position = 1
        end_position = len(contig_seqs[contig_id])
        strand = 0
        new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                            end_position,
                                                            strand),
                                                            #),
                                                            type="gene")
        new_feature_gene.qualifiers['locus_tag'] = id_gene # + "_0001"
        # print(new_feature_gene.qualifiers['locus_tag'] )
        # Add gene information to contig record.
        record.features.append(new_feature_gene)

        new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                end_position,
                                                                # strand),
                                                                ),
                                                            type="CDS")

        new_feature_cds.qualifiers['translation'] = gene_protein_seq[contig_id] #ad_gene
        new_feature_cds.qualifiers['locus_tag'] = id_gene # + "_0001"

        # Add GO annotation according to the namespace.
        # print(contig_id)
        if contig_id in annotation_data.keys():
            gene_gos = re.split(';|,', annotation_data[contig_id]['GOs'])
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
        if contig_id in annotation_data.keys():
            gene_ecs = re.split(';|,', annotation_data[contig_id]['EC'])
            if gene_ecs != [""]:
                new_feature_cds.qualifiers['EC_number'] = [ec.replace('ec:', '') for ec in gene_ecs]

        # Add CDS information to contig record
        record.features.append(new_feature_cds)
        # print(record.features)
        seq_objects.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(seq_objects, gbk_out, 'genbank')

def main(genome_fasta, prot_fasta, annot_table, species_name, gbk_out, gobasic_file = None):
    # check validity of inputs
    for elem in [genome_fasta, prot_fasta]:
        if not is_valid_file(elem):
            logger.critical(f"{elem} is not a valid path file.")
            sys.exit(1)
    if gobasic_file:
        if not is_valid_file(gobasic_file):
            logger.critical(f"{gobasic_file} is not a valid path file.")
            sys.exit(1)


    faa_to_gbk(genome_fasta, prot_fasta, annot_table, species_name, gbk_out, gobasic_file)
