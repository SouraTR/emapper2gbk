#!/usr/bin/env python

"""Tests for `eggnog2gbk` package."""

import os
import shutil
import subprocess
from eggnog2gbk.eggnog2gbk import gbk_creation
from Bio import SeqIO


GFF_INPUT = "betaox.gff"
FAA_INPUT = "betaox.faa"
FAA_DIR = "faa/"
FNA_INPUT = "betaox.fna"
FNA_DIR = "fna/"
ANNOT_INPUT = "betaox_annotation.tsv"
ANNOT_DIR = "ann/"
ORG_NAME = "Escherichia coli"
ORG_FILE = "organism_names.tsv"

EXPECTED_GBK_NO_GFF = "betaox_no_gff.gbk"
EXPECTED_GBK_WITH_GFF = "betaox_from_gff.gbk"

ANNOTATIONS_TYPES = ['go_function', 'go_process', 'go_component', 'EC_number']

ANNOTATIONS_BY_GENOME = {'gene1781':{'go_component':['GO:0005575'],
                                    'go_function':['GO:0003674'],
                                    'go_process':[],
                                    'EC_number':[]
                                    },
                        'gene1887':{'go_component':['GO:0005575'],
                                    'go_function':['GO:0003674'],
                                    'go_process':['GO:0008150'],
                                    'EC_number':['6.2.1.3']
                                    },
                        'gene2441':{'go_function':['GO:0003674'],
                                    'go_process':['GO:0008150'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    },
                        'gene3987':{'go_component':['GO:0005575'],
                                    'go_function':['GO:0003674'],
                                    'go_process':['GO:0008150'],
                                    'EC_number':['2.3.1.16']
                                    },
                        'gene3988':{'go_component':[],
                                    'go_function':['GO:0003674'],
                                    'go_process':['GO:0008150'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    }
                        }

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

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
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

    compare_two_gbks(EXPECTED_GBK_WITH_GFF, gbk_test)
    os.remove(gbk_test)
    return

def compare_two_gbks(expected_gbk:str, tested_gbk:str):
    """Compare the annotations of 2 genbank files.

    Args:
        expected_gbk (str): path to expected gbk
        tested_gbk (str): path to the second gbk
    """
    loaded_gbk = SeqIO.to_dict(SeqIO.parse(expected_gbk, "genbank"))
    loaded_test = SeqIO.to_dict(SeqIO.parse(tested_gbk, "genbank"))

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
    return

def gbk_from_dir_test():
    """Test genomic mode with directories.
    """
    gbk_dir_test = 'gbk_g/'
    os.makedirs(gbk_dir_test)
    gbk_creation(genome=FNA_DIR,
                proteome=FAA_DIR,
                annot=ANNOT_DIR,
                org=ORG_FILE,
                gff=None,
                gbk=gbk_dir_test,
                gobasic='go-basic.obo',
                dirmode=True)

    for gbk in os.listdir(gbk_dir_test):
        loaded_gbk = SeqIO.to_dict(SeqIO.parse(f"{gbk_dir_test}/{gbk}", "genbank"))
        annotations = {i:None 
                        for i in loaded_gbk.keys()}
        for gene in annotations: 
            for index in range(0, len(loaded_gbk[gene].features)):
                if loaded_gbk[gene].features[index].type=='CDS':
                    annotations[gene] = loaded_gbk[gene].features[index].qualifiers
                    # check annotations
                    for ann in ANNOTATIONS_TYPES:
                        if ann in annotations[gene]:
                            print(set(annotations[gene][ann]), set(ANNOTATIONS_BY_GENOME[gene][ann]))
                            assert set(annotations[gene][ann]) == set(ANNOTATIONS_BY_GENOME[gene][ann])

    shutil.rmtree(gbk_dir_test)
    return

def gbk_metagenomic_mode_test():
    """Test metagenomic mode.
    """
    gbk_dir_test = 'gbk_mg/'
    #TODO
    shutil.rmtree(gbk_dir_test)
    return
    
if __name__ == "__main__":
    gbk_no_gff_test()
    # gbk_from_gff_test()
    # gbk_from_dir_test()
    # gbk_metagenomic_mode_test()