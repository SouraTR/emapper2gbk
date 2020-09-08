import sys
import os
import pandas as pa
import csv
import datetime
import itertools
import logging
import numpy as np
import pronto
import requests
import shutil

from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature as sf

try:
    # Import to be compatible with biopython version lesser than 1.78
    from Bio.Alphabet import IUPAC
except ImportError:
    # Exception to be compatible with biopython version superior to 1.78
    IUPAC = None

logger = logging.getLogger(__name__)


def get_basename(filepath):
    """Return the basename of given filepath.
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: basename

    >>> basename('~/an/interesting/file.txt')
    'file
    """
    return os.path.splitext(os.path.basename(filepath))[0]


def get_extension(filepath):
    """Get the extension of a filepath
    
    Args:
        filepath (str): path to a file
    
    Returns:
        str: extention of the file

    >>> extension('~/an/interesting/file.lp')
    'lp'
    >>> extension('nothing')
    ''
    >>> extension('nothing.important')
    'important'
    """
    return os.path.splitext(os.path.basename(filepath))[1][1:]

def is_valid_path(filepath):
    """Return True if filepath is valid.
    
    Args:
        filepath (str): path to file
    
    Returns:
        bool: True if path exists, False otherwise
    """
    if filepath and not os.access(filepath, os.W_OK):
        try:
            open(filepath, 'w').close()
            os.unlink(filepath)
            return True
        except OSError:
            return False
    else:  # path is accessible
        return True


def is_valid_file(filepath):
    """Return True if filepath exists.

    Args:
        filepath (str): path to file

    Returns:
        bool: True if path exists, False otherwise
    """
    try:
        open(filepath, 'r').close()
        return True
    except OSError:
        return False


def is_valid_dir(dirpath):
    """Return True if directory exists or can be created (then create it).
    
    Args:
        dirpath (str): path of directory

    Returns:
        bool: True if dir exists, False otherwise
    """
    if not os.path.isdir(dirpath):
        try:
            os.makedirs(dirpath)
            return True
        except OSError:
            return False
    else:
        return True


def create_GO_namespaces_alternatives(gobasic_file = None):
    """
    Use pronto to query the Gene Ontology and to create the Ontology.
    Create a dictionary which contains for all GO terms their GO namespaces (molecular_function, ..).
    Create a second dictionary containing alternative ID for some GO terms (deprecated ones).
    """
    go_basic_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'

    if gobasic_file:
        if is_valid_file(gobasic_file):
            go_ontology = pronto.Ontology(gobasic_file)
        else:
            logger.critical(gobasic_file + ' file is not available, emapper2gbk will download the GO ontology and create it.')
            response = requests.get(go_basic_obo_url, stream=True)
            with open(gobasic_file, 'wb') as go_basic_file_write:
                shutil.copyfileobj(response.raw, go_basic_file_write)
            go_ontology = pronto.Ontology(go_basic_obo_url)
    else:
        go_ontology = pronto.Ontology(go_basic_obo_url)

    # For each GO terms look to the namespaces associated with them.
    go_namespaces = {}
    for go_term in go_ontology:
        if 'GO:' in go_term:
            go_namespaces[go_term] = go_ontology[go_term].namespace

    # For each GO terms look if there is an alternative ID fo them.
    go_alternative = {}
    for go_term in go_ontology:
        if go_ontology[go_term].alternate_ids != frozenset():
            for go_alt in go_ontology[go_term].alternate_ids:
                go_alternative[go_alt] = go_term

    return go_namespaces, go_alternative


def create_taxonomic_data(species_name):
    """
    Query the EBI with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.
    """
    species_informations = {}
    species_name_url = species_name.replace(' ', '%20')

    if species_name == "bacteria":
        species_informations = {'db_xref': 'taxon:2', 'scientificName': 'Bacteria', 'commonName': 'eubacteria', 'formalName': 'false', 'rank': 'superkingdom', 'data_file_division': 'PRO', 'geneticCode': '11', 'submittable': 'false', 'description': 'bacteria genome', 'organism': 'bacteria', 'keywords': ['bacteria']} 
    elif species_name == "metagenome":
        species_informations = {'db_xref': 'taxon:256318', 'scientificName': 'metagenome', 'formalName': 'false', 'rank': 'species', 'division': 'UNC', 'lineage': 'unclassified sequences; metagenomes; ', 'geneticCode': '11', 'mitochondrialGeneticCode': '2', 'plastIdGeneticCode': '11', 'submittable': 'true'}
    elif species_name == 'cellular organisms':
        species_informations = {'db_xref': 'taxon:131567', 'scientificName': 'cellular organisms', 'formalName': 'false', 'rank': 'no rank', 'division': 'UNC', 'geneticCode': '1', 'submittable': 'false'}
    else:
        url = 'https://www.ebi.ac.uk/ena/data/taxonomy/v1/taxon/scientific-name/' + species_name_url
        response = requests.get(url)
        temp_species_informations = response.json()[0]
        # print(temp_species_informations)
        for temp_species_information in temp_species_informations:
            # print(temp_species_information)
            # print(temp_species_informations[temp_species_information])
            if temp_species_information == 'lineage':
                species_informations['taxonomy'] = temp_species_informations[temp_species_information].split('; ')[:-1]
            elif temp_species_information == 'division':
                species_informations['data_file_division'] = temp_species_informations[temp_species_information]
            elif temp_species_information == 'taxId':
                species_informations['db_xref'] = 'taxon:' + str(temp_species_informations[temp_species_information])
            else:
                species_informations[temp_species_information] = temp_species_informations[temp_species_information]

    compatible_species_name = species_name.replace('/', '_')
    species_informations['description'] = compatible_species_name + ' genome'
    species_informations['organism'] = compatible_species_name
    species_informations['keywords'] = [compatible_species_name]

    return species_informations


def record_info(record_id, record_seq, species_informations):
    """
    Create SeqBio record informations from species_informations dictionary and record id and record seq.
    """
    if record_id.isnumeric():
        newname = f"_{record_id}"
    elif "|" in record_id:
        newname = record_id.split("|")[0]
    else:
        newname = record_id

    record = SeqRecord(record_seq, id=record_id, name=newname,
                    description=species_informations['description'],
                    annotations={"molecule_type": "DNA"})

    # Condition to be compatible with biopython version lesser than 1.78
    if IUPAC:
        record.seq.alphabet = IUPAC.ambiguous_dna

    if 'data_file_division' in species_informations:
        record.annotations['data_file_division'] = species_informations['data_file_division']
    record.annotations['date'] = datetime.date.today().strftime('%d-%b-%Y').upper()
    if 'topology' in species_informations:
        record.annotations['topology'] = species_informations['topology']
    record.annotations['accessions'] = record_id
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
                                                        len(record_seq)),
                                                        type="source")
    new_feature_source.qualifiers['scaffold'] = record_id
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


def read_annotation(eggnog_outfile:str):
    """Read an eggno-mapper annotation file and retrieve EC numbers and GO terms by genes.

    Args:
        eggnog_outfile (str): path to eggnog-mapper annotation file

    Returns:
        dict: dict of genes and their annotations as {gene1:{EC:'..,..', GOs:'..,..,'}}
    """
    # Retrieve headers name at line 4.
    colnames_linenb = 3
    with open(eggnog_outfile, 'r') as f:
        headers_row = next(itertools.islice(csv.reader(f), colnames_linenb, None))[0].lstrip("#").strip().split('\t')

    # Fix issue when header is incomplete (eggnog before version 2.0).
    if len(headers_row) == 17:
        headers_row.extend(['tax_scope', 'eggNOG_OGs', 'bestOG', 'COG_functional_category', 'eggNOG_free_text'])

    # Use chunk when reading eggnog file to cope with big file.
    chunksize = 10 ** 6
    for annotation_data in pa.read_csv(eggnog_outfile, sep='\t', comment='#', header=None, dtype = str, chunksize = chunksize):
        annotation_data.replace(np.nan, '', inplace=True)
        # Assign the headers
        annotation_data.columns = headers_row
        annotation_dict = annotation_data.set_index('query_name')[['GOs','EC', 'Preferred_name']].to_dict('index')
        for key in annotation_dict:
            yield key, annotation_dict[key]


def create_cds_feature(id_gene, start_position, end_position, strand, annotation_data, go_namespaces, go_alternatives, gene_protein_seq):
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

    return new_feature_cds
