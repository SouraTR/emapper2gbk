import sys
import os

def is_valid_path(filepath):
    """Return True if filepath is valid
    
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
    """Return True if filepath exists

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
    """Return True if directory exists or can be created (then create it)
    
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
            go_namespaces[go_term] = go_ontology[go_term].name
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
