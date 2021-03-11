# Copyright (C) 2019-2021 Clémence Frioux & Arnaud Belcour - Inria Dyliss - Pleiade
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import sys
import argparse
import os
import re
import shutil
import logging

from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
from typing import Union

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


def faa_to_gbk(nucleic_fasta:str, protein_fasta:str, annotation_data:Union[str, dict], species_name:str, gbk_out:str, gobasic:Union[None, str, dict], merge:int):
    """
    From a genome fasta (containing each genes of the genome),
    a protein fasta (containing each protein sequence),
    an eggnog annotation table,
    a species name,
    a name for the genbank output,
    gobasic: name of go-basic.obo or a tuple of dictionaries (go_namespaces, go_alternatives),
    create a genbank file.
    """
    genome_id = get_basename(nucleic_fasta)

    logger.info('Formatting fasta and annotation file for ' + genome_id)

    # Dictionary with gene id as key and nucleic sequence as value.
    gene_nucleic_seqs = OrderedDict()

    for record in SeqIO.parse(nucleic_fasta, "fasta"):
        gene_nucleic_seqs[record.id] = record.seq

    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seqs = OrderedDict()

    for record in SeqIO.parse(protein_fasta, "fasta"):
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

    logger.info('Assembling Genbank informations for ' + genome_id)

    # Create fake contig by merging genes.
    if merge:
        create_genbank_fake_contig(gene_nucleic_seqs, gene_protein_seqs, annotation_data, go_namespaces, go_alternatives, gbk_out, species_informations, merge)
    else:
        create_genbank(gene_nucleic_seqs, gene_protein_seqs, annotation_data, go_namespaces, go_alternatives, gbk_out, species_informations)


def create_genbank(gene_nucleic_seqs, gene_protein_seqs, annotation_data, go_namespaces, go_alternatives, gbk_out, species_informations):
    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    records = []

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
        new_feature_cds.qualifiers['locus_tag'] = id_gene

        # Add functional annotations.
        if gene_nucleic_id in annotation_data.keys():
            # Add gene name.
            if 'Preferred_name' in annotation_data[gene_nucleic_id]:
                if annotation_data[gene_nucleic_id]['Preferred_name'] != "":
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
                if '' in gene_ecs:
                    gene_ecs.remove('')
                if '-' in gene_ecs:
                    gene_ecs.remove('-')
                if gene_ecs != []:
                    new_feature_cds.qualifiers['EC_number'] = gene_ecs

        if gene_nucleic_id in gene_protein_seqs.keys():
            # Add protein sequence.
            new_feature_cds.qualifiers['translation'] = gene_protein_seqs[gene_nucleic_id]

        # Add CDS information to contig record
        record.features.append(new_feature_cds)

        records.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(records, gbk_out, 'genbank')


def create_genbank_fake_contig(gene_nucleic_seqs, gene_protein_seqs, annotation_data, go_namespaces, go_alternatives, gbk_out, species_informations, gene_per_fake_contig):
    fake_contigs = {}
    # Iterate through each contig/gene.
    maximal_contig = int(len(list(gene_nucleic_seqs.keys())) / gene_per_fake_contig)

    maximal_contig_str_stize = len(str(maximal_contig))

    def homogonize_id(contig_id, max_id_length):
        old_contig_number= contig_id.split('_')[-1]
        if len(old_contig_number) < max_id_length:
            nb_diff = max_id_length - len(old_contig_number)
            new_contig_id = '_'.join(contig_id.split('_')[:-1]) + '_' + str(nb_diff*[0]) + str(old_contig_number)
        else:
            new_contig_id = contig_id

        return new_contig_id

    fake_contigs.update({homogonize_id('fake_contig_'+str(index), maximal_contig_str_stize):list(gene_nucleic_seqs.keys())[x:x+gene_per_fake_contig] for index, x in enumerate(range(0, len(gene_nucleic_seqs), gene_per_fake_contig))})

    # All SeqRecord objects will be stored in a list and then give to the SeqIO writer to create the genbank.
    records = []
    for contig in fake_contigs:
        contig_seq = Seq('').join([gene_nucleic_seqs[gene] for gene in fake_contigs[contig]])

        # Create a SeqRecord object using gene information.
        record = record_info(contig, contig_seq, species_informations)

        gene_position = 0
        for id_gene in fake_contigs[contig]:
            start_position = gene_position
            gene_position += len(gene_nucleic_seqs[id_gene])
            end_position = gene_position
            strand = 0
            new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                end_position,
                                                                strand),
                                                                type="gene")
            new_feature_gene.qualifiers['locus_tag'] = id_gene

            # Add gene information to contig record.
            record.features.append(new_feature_gene)

            new_feature_cds = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                end_position,
                                                                strand),
                                                                type="CDS")

            # Add gene ID in locus_tag.
            new_feature_cds.qualifiers['locus_tag'] = id_gene

            # Add functional annotations.
            if id_gene in annotation_data.keys():
                # Add gene name.
                if 'Preferred_name' in annotation_data[id_gene]:
                    if annotation_data[id_gene]['Preferred_name'] != "":
                        new_feature_cds.qualifiers['gene'] = annotation_data[id_gene]['Preferred_name']

                # Add GO annotation according to the namespace.
                if 'GOs' in annotation_data[id_gene]:
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
                    if '' in gene_ecs:
                        gene_ecs.remove('')
                    if '-' in gene_ecs:
                        gene_ecs.remove('-')
                    if gene_ecs != []:
                        new_feature_cds.qualifiers['EC_number'] = gene_ecs

            if id_gene in gene_protein_seqs.keys():
                # Add protein sequence.
                new_feature_cds.qualifiers['translation'] = gene_protein_seqs[id_gene]

            # Add CDS information to contig record
            record.features.append(new_feature_cds)

        records.append(record)
    SeqIO.write(records, gbk_out, 'genbank')

def main(nucleic_fasta, protein_fasta, annotation_data, species_name, gbk_out, gobasic=None, merge=None):
    # check validity of inputs
    for elem in [nucleic_fasta, protein_fasta]:
        if not is_valid_file(elem):
            logger.critical(f"{elem} is not a valid path file.")
            sys.exit(1)

    faa_to_gbk(nucleic_fasta, protein_fasta, annotation_data, species_name, gbk_out, gobasic, merge)
