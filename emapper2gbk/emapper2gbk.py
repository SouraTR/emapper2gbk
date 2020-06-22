"""Main module."""
import logging
import os
import csv
import sys
from multiprocessing import Pool
from emapper2gbk import gff_to_gbk
from emapper2gbk import fa_to_gbk
from emapper2gbk.utils import create_GO_namespaces_alternatives, get_basename, is_valid_file, get_extension, read_annotation

logger = logging.getLogger(__name__)

def gbk_creation(genome:str, proteome:str, annot:str, org:str, gbk:str, gobasic:str, dirmode:bool=False, cpu:int=1, metagenomic_mode:bool=False, gff:str=None):
    """Create gbk files from genomic information and eggnog-mapper annotation outputs.

    Args:
        genome (str): genome file or dir
        proteome (str): proteome file or dir
        annot (str): annotation file or dir
        org (str): organims name or mapping file
        gbk (str): output file or directory
        gobasic (str): path to go-basic.obo file
        dirmode (bool, optional): directory mode (instead of single file). Defaults to False.
        dirmode (bool, optional): metagenomic mode. Defaults to False.
        cpu (int, optional): number of cpu, used for multi process in directory mode. Defaults to 1.
        gff (str, optional): gff file or dir. Defaults to None.
    """
    if not dirmode:
        if gff:
            gff_to_gbk.main(genome, proteome, annot, gff, org, gbk, gobasic)
        else:
            fa_to_gbk.main(genome, proteome, annot, org, gbk, gobasic)
    else:
        gbk_pool = Pool(processes=cpu)

        all_genomes = set([get_basename(i) for i in os.listdir(genome)])
        # filename or single org name
        if is_valid_file(org):
            with open(org, 'r') as csvfile:
                try:
                    dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters=";,\t")
                except csv.Error:
                    logger.critical(f"Could not determine the delimiter in the organism tabulated file.")
                    exit(1)
                csvfile.seek(0)
                reader = csv.reader(csvfile, dialect)
                org_mapping = {i[0]:i[1] for i in reader}
        else:
            org_mapping = {i:org for i in all_genomes}

        # check that all data is here
        try:
            assert set(org_mapping.keys()) == all_genomes
        except AssertionError:
            logger.critical(f"Genomes in {genome} do not match the genomes IDs of {org} (first column, no extension to the genome names).")
            sys.exit(1)
        try:
            assert all(ext == 'fna' for ext in [get_extension(i) for i in os.listdir(genome)])
        except AssertionError:
            logger.critical(f"Genomes dir {genome} must contain only '.fna' files")
            sys.exit(1)
        #TODO later: be less hard on that constraint, just check that every genome foo.fna has a foo.faa in proteome, a foo.tsv in annot etc. This could make possible for a user to give the same dir for all inputs. Or ignore the extension and ensure thee is a foo.* in annot, a foo.* in proteome etc.
        try:
            set([get_basename(i) for i in os.listdir(proteome)]) == all_genomes
        except AssertionError:
            logger.critical(f"Genomes names in {genome} do not match the names in {proteome}.")
            sys.exit(1)
        try:
            assert all(ext == 'faa' for ext in [get_extension(i) for i in os.listdir(proteome)])
        except AssertionError:
            logger.critical(f"Proteomes dir {proteome} must contain only '.faa' files")
            sys.exit(1)
        if os.path.isdir(annot):
            try:
                assert set([get_basename(i) for i in os.listdir(annot)]) == all_genomes
            except AssertionError:
                logger.critical(f"Genomes names in {genome} do not match the names in {annot}.")
                sys.exit(1)
            try:
                assert all(ext == 'tsv' for ext in [get_extension(i) for i in os.listdir(annot)])
            except AssertionError:
                logger.critical(f"Annotations dir {annot} must contain only '.tsv' files")
                sys.exit(1)
        if gff:
            try:
                assert set([get_basename(i) for i in os.listdir(gff)]) == all_genomes
            except AssertionError:
                logger.critical(f"Genomes names in {genome} do not match the names in {gff}.")
                sys.exit(1)
            try:
                assert all(ext == 'gff' for ext in [get_extension(i) for i in os.listdir(gff)])
            except AssertionError:
                logger.critical(f"GFF dir {gff} must contain only '.gff' files")
                sys.exit(1)

        multiprocess_data = []

        # Query Gene Ontology to extract namespaces and alternative IDs.
        # go_namespaces: Dictionary GO id as term and GO namespace as value.
        # go_alternatives: Dictionary GO id as term and GO alternatives id as value.
        go_namespaces, go_alternatives = create_GO_namespaces_alternatives(gobasic)

        if gff and not metagenomic_mode:
            for genome_id in all_genomes:
                    multiprocess_data.append(
                        {'genome':f"{genome}/{genome_id}.fna",
                        'proteome':f"{proteome}/{genome_id}.faa",
                        'annot':f"{annot}/{genome_id}.tsv",
                        'gff':f"{gff}/{genome_id}.gff",
                        'org':org_mapping[genome_id],
                        'gbk':f"{gbk}/{genome_id}.gbk",
                        'gobasic':(go_namespaces, go_alternatives)}
                        )
            gbk_pool.map(run_gff_to_gbk, multiprocess_data)
        elif gff and metagenomic_mode:
            # read annotation of gene catalogue
            annot_genecat = dict(read_annotation(annot))
            for genome_id in all_genomes:
                    multiprocess_data.append(
                        {'genome':f"{genome}/{genome_id}.fna",
                        'proteome':f"{proteome}/{genome_id}.faa",
                        'annot':annot_genecat,
                        'gff':f"{gff}/{genome_id}.gff",
                        'org':org_mapping[genome_id],
                        'gbk':f"{gbk}/{genome_id}.gbk",
                        'gobasic':(go_namespaces, go_alternatives)}
                        )
            gbk_pool.map(run_gff_to_gbk, multiprocess_data)
        elif not gff and not metagenomic_mode:
            for genome_id in all_genomes:
                    multiprocess_data.append(
                        {'genome':f"{genome}/{genome_id}.fna",
                        'proteome':f"{proteome}/{genome_id}.faa",
                        'annot':f"{annot}/{genome_id}.tsv",
                        'org':org_mapping[genome_id],
                        'gbk':f"{gbk}/{genome_id}.gbk",
                        'gobasic':(go_namespaces, go_alternatives)}
                        )
            gbk_pool.map(run_fa_to_gbk, multiprocess_data)
        else:
            # read annotation of gene catalogue
            annot_genecat = dict(read_annotation(annot))
            for genome_id in all_genomes:
                    multiprocess_data.append(
                        {'genome':f"{genome}/{genome_id}.fna",
                        'proteome':f"{proteome}/{genome_id}.faa",
                        'annot':annot_genecat,
                        'org':org_mapping[genome_id],
                        'gbk':f"{gbk}/{genome_id}.gbk",
                        'gobasic':(go_namespaces, go_alternatives)}
                        )
            gbk_pool.map(run_fa_to_gbk, multiprocess_data)

        gbk_pool.close()
        gbk_pool.join()


def run_gff_to_gbk(mp_data):
    """Run gff_to_gbk from a dict of arguments.

    Args:
        mp_data (dict): all arguments
    """
    logger.info(get_basename(mp_data["genome"]))
    gff_to_gbk.main(mp_data["genome"], 
                    mp_data["proteome"], 
                    mp_data["annot"], 
                    mp_data["gff"], 
                    mp_data["org"], 
                    mp_data["gbk"], 
                    mp_data["gobasic"]
                    )


def run_fa_to_gbk(mp_data):
    """Run fa_to_gbk from a dict of arguments.

    Args:
        mp_data (dict): all arguments
    """
    logger.info(get_basename(mp_data["genome"]))
    fa_to_gbk.main(mp_data["genome"], 
                    mp_data["proteome"], 
                    mp_data["annot"], 
                    mp_data["org"], 
                    mp_data["gbk"], 
                    mp_data["gobasic"]
                    )