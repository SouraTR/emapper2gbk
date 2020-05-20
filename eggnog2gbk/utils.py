import sys
import os
import pandas as pa
import csv
import itertools
import numpy as np
import pronto
import requests

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

def create_GO_dataframes(gobasic_file = None):
    """
    Use pronto to query the Gene Ontology and to create the Ontology.
    Create a dataframe which contains for all GO terms their GO namespaces (molecular_function, ..).
    Create a second dataframe containing alternative ID for some GO terms (deprecated ones).
    """
    if gobasic_file:
        go_ontology = pronto.Ontology(gobasic_file)
    else:
        go_ontology = pronto.Ontology('http://purl.obolibrary.org/obo/go/go-basic.obo')

    # For each GO terms look to the namespaces associated with them.
    go_namespaces = {}
    for go_term in go_ontology:
        if 'GO:' in go_term:
            go_namespaces[go_term] = go_ontology[go_term].namespace
    df_go_namespace = pa.DataFrame.from_dict(go_namespaces, orient='index')
    df_go_namespace.reset_index(inplace=True)
    df_go_namespace.columns = ['GO', 'namespace']

    # For each GO terms look if there is an alternative ID fo them.
    go_alt_ids = {}
    for go_term in go_ontology:
        if go_ontology[go_term].alternate_ids != frozenset():
            for go_alt in go_ontology[go_term].alternate_ids:
                go_alt_ids[go_alt] = go_term
    df_go_alternative = pa.DataFrame.from_dict(go_alt_ids, orient='index')
    df_go_alternative.reset_index(inplace=True)
    df_go_alternative.columns = ['GO', 'alternative_GO']

    return df_go_namespace, df_go_alternative

def create_taxonomic_data(species_name):
    """
    Query the EBI with the species name to create a dictionary containing taxon id,
    taxonomy and some other informations.
    """
    species_informations = {}
    species_name_url = species_name.replace(' ', '%20')

    if species_name == "bacteria":
        species_informations = {'db_xref': 'taxon:2', 'scientificName': 'Bacteria', 'commonName': 'eubacteria', 'formalName': 'false', 'rank': 'superkingdom', 'data_file_division': 'PRO', 'geneticCode': '11', 'submittable': 'false', 'description': 'bacteria genome', 'organism': 'bacteria', 'keywords': ['bacteria']} 
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

def read_annotation(eggnog_outfile:str):
    """Read an eggno-mapper annotation file and retrieve EC numbers and GO terms by genes.

    Args:
        eggnog_outfile (str): path to eggnog-mapper annotation file

    Returns:
        dict: dict of genes and their annotations as {gene1:{EC:'..,..', GOs:'..,..,'}}
    """
    # read eggnog-mapper annotation file
    # warning1: column names are commented
    # warning2: there are less column names than actual columns
    annotation_data = pa.read_csv(eggnog_outfile, sep='\t', comment='#', header=None, dtype = str)
    annotation_data.replace(np.nan, '', inplace=True)
    # retrieve headers name at line 4
    colnames_linenb = 3
    with open(eggnog_outfile, 'r') as f:
        headers_row = next(itertools.islice(csv.reader(f), colnames_linenb, None))[0].lstrip("#").strip().split('\t')
    # fill headers with unknowns up to columns number
    i = 1
    while len(headers_row) < len(annotation_data.columns):
        headers_row.append(f"unknown_{i}")
        i += 1

    # assign the headers:
    annotation_data.columns = headers_row

    # now create the dict
    emapper_dict = annotation_data.set_index('query_name')[['GOs','EC']].to_dict('index') # annotation_data.groupby("query_name")[["GOs","EC"]].apply(lambda x: dict(x.values)).to_dict()

    return emapper_dict
