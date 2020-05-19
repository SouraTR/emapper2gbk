#!/usr/bin/env python

"""Tests for `eggnog2gbk` package."""

import os
import shutil
import subprocess
from eggnog2gbk.eggnog2gbk import gbk_creation
from Bio import SeqIO


GFF_INPUT = "betaox.gff"
FAA_INPUT = "betaox.faa"
FNA_INPUT = "betaox.fna"
ANNOT_INPUT = "betaox_annotation.tsv"
ORG_NAME = "Escherichia coli"

EXPECTED_GBK = "betaox_no_gff.gbk"

ANNOTATIONS_TYPES = ['go_function', 'go_process', 'EC_number']

def gbk_no_gff_test():
    """Test genomic mode without gff as input.
    """
    gbk_test = 'test_no_gff.gbk'
    gbk_creation(genome=FNA_INPUT,
                proteome=FAA_INPUT,
                annot=ANNOT_INPUT,
                org=ORG_NAME,
                gbk=gbk_test,
                gff=None,
                gobasic='go-basic.obo')

    loaded_gbk = SeqIO.to_dict(SeqIO.parse(EXPECTED_GBK, "genbank"))
    loaded_test = SeqIO.to_dict(SeqIO.parse(gbk_test, "genbank"))

    assert set(loaded_gbk.keys()) == set(loaded_test.keys())

    # annotations_expected = {i:loaded_gbk[i].features[e].qualifiers 
    #                             if loaded_gbk[i].features[e].type=='CDS' 
    #                             else None 
    #                             for e in range(0, len(loaded_gbk[i].features)) 
    #                         for i in loaded_gbk.keys()}
    
    # annotations_test = {i:loaded_test[i].features[e].qualifiers 
    #                             if loaded_test[i].features[e].type=='CDS' 
    #                             else None 
    #                             for e in range(0, len(loaded_test[i].features)) 
    #                         for i in loaded_test.keys()}

    annotations_expected = {i:None 
                            for i in loaded_gbk.keys()}
    
    annotations_test = {i:None 
                        for i in loaded_test.keys()}

    for gene in annotations_expected: 
        for index in range(0, len(loaded_gbk[gene].features)):
            if loaded_gbk[gene].features[index].type=='CDS':
                annotations_expected[gene] = loaded_gbk[gene].features[index].qualifiers

    for gene in annotations_test: 
        for index in range(0, len(loaded_test[gene].features)):
            if loaded_test[gene].features[index].type=='CDS':
                annotations_test[gene] = loaded_test[gene].features[index].qualifiers

    for gene in annotations_expected:
        for qualifier in annotations_expected[gene]:
            if qualifier in ANNOTATIONS_TYPES:
                assert set(annotations_expected[gene][qualifier]) == set(annotations_test[gene][qualifier])

    os.remove(gbk_test)
    return


def gbk_from_gff_test():
    """Test genomic mode with a gff as input.
    """
    gbk_test = 'test_gff.gbk'

    gbk_creation(genome=FNA_INPUT,
                proteome=FAA_INPUT,
                annot=ANNOT_INPUT,
                gff=GFF_INPUT,
                org=ORG_NAME,
                gbk=gbk_test,
                gobasic='go-basic.obo')

    loaded_gbk = SeqIO.to_dict(SeqIO.parse(EXPECTED_GBK, "genbank"))
    loaded_test = SeqIO.to_dict(SeqIO.parse(gbk_test, "genbank"))

    os.remove(gbk_test)
    #TODO
    return

def gbk_from_dir_test():
    """Test genomic mode with directories.
    """
    #TODO
    return

def gbk_metagenomic_mode_test():
    """Test metagenomic mode.
    """
    #TODO
    return
    
if __name__ == "__main__":
    gbk_no_gff_test()
    # gbk_from_gff_test()
    # gbk_from_dir_test()
    # gbk_metagenomic_mode_test()