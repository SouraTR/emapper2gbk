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

import gffutils
import logging
import re
import sys

from Bio import SeqFeature as sf
from Bio import SeqIO
from collections import OrderedDict
from emapper2gbk.utils import check_valid_path, is_valid_file, create_GO_namespaces_alternatives, read_annotation, create_taxonomic_data, get_basename, record_info, create_cds_feature
from typing import Union

logger = logging.getLogger(__name__)


def strand_change(input_strand):
    """
    The input is strand in str ('-', '+') modify it to be a strand in int (-1, +1) to 
    be compatible with SeqIO strand reading.

    Args:
        input_strand (str): input strand
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


def gff_to_gbk(nucleic_fasta:str, protein_fasta:str, annot:Union[str, dict],
                gff:str, gff_type:str, org:str, output_path:str, gobasic:Union[None, str, dict],
                keep_gff_annot:Union[None,bool]):
    """ Create genbank file from nucleic, protein fasta and a gff file plus eggnog mapper annotation file.

    Args:
        nucleic_fasta (str): nucleic fasta file
        protein_fasta (str): protein fasta file
        annot (str): annotation file or dictionary
        gff (str): gff file
        gff_type (str): format of the gff file (default, prodigal)
        org (str): organims name or mapping file
        output_path (str): output file or directory
        gobasic (str, dict): path to go-basic.obo file or dictionary
        keep_gff_annot (bool): copy the annotation present in the GFF file into the Genbank file.
    """
    check_valid_path([nucleic_fasta, protein_fasta, gff])

    genome_id = get_basename(nucleic_fasta)

    logger.info('Creating GFF database (gffutils) for ' + genome_id)
    # Create the gff database file.
    # gffutils use sqlite3 file-based database to access data inside GFF.
    # ':memory:' ask gffutils to keep database in memory instead of writting in a file.
    gff_database = gffutils.create_db(gff, ':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

    logger.info('Formatting fasta and annotation file for ' + genome_id)
    # Dictionary with region id (contig, chromosome) as key and sequence as value.
    genome_nucleic_sequence = OrderedDict()
    for record in SeqIO.parse(nucleic_fasta, "fasta"):
        region_id = record.id
        genome_nucleic_sequence[region_id] = record.seq

    # Dictionary with gene id as key and protein sequence as value.
    gene_protein_seqs = {}

    for record in SeqIO.parse(protein_fasta, "fasta"):
        gene_protein_seqs[record.id] = record.seq

    # Create a taxonomy dictionary querying the EBI.
    species_informations = create_taxonomic_data(org)

    # Read the eggnog tsv file containing GO terms and EC associated with gene name.
    # if metagenomic mode, annotation is already read and given as a dict
    if not type(annot) is dict:
        annot = dict(read_annotation(annot))

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

    annotations_in_gff = ['product']
    # Iterate through each contig.
    # Then iterate through gene and throug RNA linked with the gene.
    # Then look if protein informations are available.
    for region_id in genome_nucleic_sequence:
        record = record_info(region_id, genome_nucleic_sequence[region_id], species_informations)
        if gff_type == 'default':
            gene_region_id = [gene for gene in gff_database.features_of_type('gene') if gene.chrom == region_id]
            for gene in gene_region_id:
                id_gene = gene.id

                # If id is numeric, change it
                if id_gene.isnumeric():
                    id_gene = f"gene_{id_gene}"
                else:
                    id_gene = id_gene

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
                    start_position = cds_object.start - 1
                    end_position = cds_object.end
                    strand = strand_change(cds_object.strand)

                    if keep_gff_annot:
                        gff_extracted_annotations = {annotation: cds_object.attributes[annotation]
                                                        for annotation in annotations_in_gff
                                                        if annotation in cds_object.attributes}
                    else:
                        gff_extracted_annotations = None

                    new_cds_feature = create_cds_feature(cds_id, start_position, end_position,
                                                        strand, annot, go_namespaces, go_alternatives,
                                                        gene_protein_seqs, gff_extracted_annotations)
                    new_cds_feature.qualifiers['locus_tag'] = cds_id
                    # Add CDS information to contig record
                    record.features.append(new_cds_feature)

        elif gff_type == 'cds_only':
            cds_region_id = [gene for gene in gff_database.features_of_type('CDS') if gene.chrom == region_id]
            for cds in cds_region_id:
                id_cds = cds.id

                # If id is numeric, change it
                if id_cds.isnumeric():
                    id_cds = f"gene_{id_cds}"
                else:
                    id_cds = id_cds

                start_position = cds.start -1
                end_position = cds.end
                strand = strand_change(cds.strand)
                new_feature_gene = sf.SeqFeature(sf.FeatureLocation(start_position,
                                                                    end_position,
                                                                    strand),
                                                                    type="gene")
                new_feature_gene.qualifiers['locus_tag'] = id_cds
                # Add gene information to contig record.
                record.features.append(new_feature_gene)

                new_cds_feature = create_cds_feature(id_cds, start_position, end_position, strand, annot, go_namespaces, go_alternatives, gene_protein_seqs)
                new_cds_feature.qualifiers['locus_tag'] = id_cds
                # Add CDS information to contig record
                record.features.append(new_cds_feature)

        seq_objects.append(record)

    # Create Genbank with the list of SeqRecord.
    SeqIO.write(seq_objects, output_path, 'genbank')
