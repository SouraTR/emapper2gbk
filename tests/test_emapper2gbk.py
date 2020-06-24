#!/usr/bin/env python

"""Tests for `emapper2gbk` package."""

import os
import shutil
import subprocess
from emapper2gbk.emapper2gbk import gbk_creation
from Bio import SeqIO


GFF_INPUT = "betaox.gff"
FAA_INPUT = "betaox.faa"
FAA_DIR = "faa/"
FNA_INPUT = "betaox.fna"
FNA_DIR = "fna/"
ANNOT_INPUT = "betaox_annotation.tsv"
ANNOT_DIR = "ann/"
GFF_DIR = "gff/"
ORG_NAME = "Escherichia coli"
ORG_FILE = "organism_names.tsv"

EXPECTED_GBK_NO_GFF = "betaox_no_gff.gbk"
EXPECTED_GBK_WITH_GFF = "betaox_from_gff.gbk"

ANNOTATIONS_TYPES = ['go_function', 'go_process', 'go_component', 'EC_number']

ANNOTATIONS_BY_GENOME = {'gene1781':{'go_component':['GO:0005575', 'GO:0005623', 'GO:0005886', 
                                                    'GO:0016020', 'GO:0044464', 'GO:0071944'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0015645', 
                                                    'GO:0016405', 'GO:0016874', 'GO:0016877', 
                                                    'GO:0016878', 'GO:0031956'],
                                    'go_process':[],
                                    'EC_number':[]
                                    },
                        'gene1887':{'go_component':['GO:0005575', 'GO:0005622', 'GO:0005623',
                                                    'GO:0005737', 'GO:0005829', 'GO:0005886',
                                                    'GO:0009898', 'GO:0016020', 'GO:0044424',
                                                    'GO:0044425', 'GO:0044444', 'GO:0044459', 
                                                    'GO:0044464', 'GO:0071944', 'GO:0098552', 
                                                    'GO:0098562'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0004467', 
                                                    'GO:0005488', 'GO:0005504', 'GO:0008289', 
                                                    'GO:0015645', 'GO:0016405', 'GO:0016874', 
                                                    'GO:0016877', 'GO:0016878', 'GO:0031406', 
                                                    'GO:0033293', 'GO:0036041', 'GO:0036094', 
                                                    'GO:0043167', 'GO:0043168', 'GO:0043177', 
                                                    'GO:0070538'],
                                    'go_process':['GO:0001676', 'GO:0006082', 'GO:0006139', 
                                                'GO:0006163', 'GO:0006629', 'GO:0006631', 
                                                'GO:0006635', 'GO:0006637', 'GO:0006644', 
                                                'GO:0006725', 'GO:0006732', 'GO:0006753', 
                                                'GO:0006790', 'GO:0006793', 'GO:0006796', 
                                                'GO:0006807', 'GO:0008150', 'GO:0008152', 
                                                'GO:0008610', 'GO:0008654', 'GO:0009056', 
                                                'GO:0009058', 'GO:0009062', 'GO:0009117', 
                                                'GO:0009150', 'GO:0009259', 'GO:0009314', 
                                                'GO:0009411', 'GO:0009416', 'GO:0009628', 
                                                'GO:0009987', 'GO:0016042', 'GO:0016054', 'GO:0019395', 'GO:0019637', 'GO:0019693', 
                                                'GO:0019752', 'GO:0030258', 'GO:0032787', 
                                                'GO:0033865', 'GO:0033875', 'GO:0034032', 
                                                'GO:0034440', 'GO:0034641', 'GO:0035383',
                                                'GO:0043436', 'GO:0043603', 'GO:0044237', 
                                                'GO:0044238', 'GO:0044242', 'GO:0044248', 
                                                'GO:0044249', 'GO:0044255', 'GO:0044281', 
                                                'GO:0044282', 'GO:0046395', 'GO:0046483', 
                                                'GO:0050896', 'GO:0051186', 'GO:0055086', 
                                                'GO:0055114', 'GO:0071704', 'GO:0072329', 
                                                'GO:0072521', 'GO:0090407', 'GO:1901135', 
                                                'GO:1901360', 'GO:1901564', 'GO:1901575', 
                                                'GO:1901576'],
                                    'EC_number':['6.2.1.3']
                                    },
                        'gene2441':{'go_function':['GO:0003674', 'GO:0003824', 'GO:0003857', 
                                                    'GO:0004165', 'GO:0004300', 'GO:0008691', 
                                                    'GO:0008692', 'GO:0016491', 'GO:0016614', 
                                                    'GO:0016616', 'GO:0016829', 'GO:0016835', 
                                                    'GO:0016836', 'GO:0016853', 'GO:0016854', 
                                                    'GO:0016856', 'GO:0016860', 'GO:0016863'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631', 
                                                    'GO:0006635', 'GO:0006725', 'GO:0006805', 
                                                    'GO:0008150', 'GO:0008152', 'GO:0009056', 
                                                    'GO:0009062', 'GO:0009404', 'GO:0009407', 
                                                    'GO:0009410', 'GO:0009636', 'GO:0009850', 
                                                    'GO:0009852', 'GO:0009987', 'GO:0010124', 
                                                    'GO:0010817', 'GO:0016042', 'GO:0016054', 
                                                    'GO:0019395', 'GO:0019439', 'GO:0019748', 
                                                    'GO:0019752', 'GO:0030258', 'GO:0032787', 
                                                    'GO:0034440', 'GO:0042178', 'GO:0042221', 
                                                    'GO:0042445', 'GO:0042447', 'GO:0042537', 
                                                    'GO:0043436', 'GO:0044237', 'GO:0044238', 
                                                    'GO:0044242', 'GO:0044248', 'GO:0044255', 
                                                    'GO:0044281', 'GO:0044282', 'GO:0046395', 
                                                    'GO:0050896', 'GO:0051716', 'GO:0055114', 
                                                    'GO:0065007', 'GO:0065008', 'GO:0070887', 
                                                    'GO:0071466', 'GO:0071704', 'GO:0072329', 
                                                    'GO:0098754', 'GO:1901360', 'GO:1901361', 
                                                    'GO:1901575'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    },
                        'gene3987':{'go_component':['GO:0005575', 'GO:0005622', 'GO:0005623', 
                                                    'GO:0005737', 'GO:0044424', 'GO:0044464'],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0003988', 
                                                    'GO:0016408', 'GO:0016740', 'GO:0016746', 
                                                    'GO:0016747'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631', 
                                                'GO:0006635', 'GO:0008150', 'GO:0008152', 
                                                'GO:0009056', 'GO:0009062', 'GO:0009987', 
                                                'GO:0016042', 'GO:0016054', 'GO:0019395', 
                                                'GO:0019752', 'GO:0030258', 'GO:0032787', 
                                                'GO:0034440', 'GO:0043436', 'GO:0044237', 
                                                'GO:0044238', 'GO:0044242', 'GO:0044248', 
                                                'GO:0044255', 'GO:0044281', 'GO:0044282', 
                                                'GO:0046395', 'GO:0055114', 'GO:0071704', 
                                                'GO:0072329', 'GO:1901575'],
                                    'EC_number':['2.3.1.16']
                                    },
                        'gene3988':{'go_component':[],
                                    'go_function':['GO:0003674', 'GO:0003824', 'GO:0003857', 
                                                'GO:0004165', 'GO:0004300', 'GO:0008691', 
                                                'GO:0008692', 'GO:0016491', 'GO:0016614', 
                                                'GO:0016616', 'GO:0016829', 'GO:0016835', 
                                                'GO:0016836', 'GO:0016853', 'GO:0016854', 
                                                'GO:0016856', 'GO:0016860', 'GO:0016863'],
                                    'go_process':['GO:0006082', 'GO:0006629', 'GO:0006631', 
                                                'GO:0006635', 'GO:0006725', 'GO:0006805', 
                                                'GO:0008150', 'GO:0008152', 'GO:0009056', 
                                                'GO:0009062', 'GO:0009404', 'GO:0009407', 
                                                'GO:0009410', 'GO:0009636', 'GO:0009850', 
                                                'GO:0009852', 'GO:0009987', 'GO:0010124', 
                                                'GO:0010817', 'GO:0016042', 'GO:0016054', 
                                                'GO:0019395', 'GO:0019439', 'GO:0019748', 
                                                'GO:0019752', 'GO:0030258', 'GO:0032787', 
                                                'GO:0034440', 'GO:0042178', 'GO:0042221', 
                                                'GO:0042445', 'GO:0042447', 'GO:0042537', 
                                                'GO:0043436', 'GO:0044237', 'GO:0044238', 
                                                'GO:0044242', 'GO:0044248', 'GO:0044255', 
                                                'GO:0044281', 'GO:0044282', 'GO:0046395', 
                                                'GO:0050896', 'GO:0051716', 'GO:0055114', 
                                                'GO:0065007', 'GO:0065008', 'GO:0070887', 
                                                'GO:0071466', 'GO:0071704', 'GO:0072329', 
                                                'GO:0098754', 'GO:1901360', 'GO:1901361', 
                                                'GO:1901575'],
                                    'EC_number':['1.1.1.157',
                                                '1.1.1.35',
                                                '4.2.1.17',
                                                '5.1.2.3',
                                                '5.3.3.8'
                                                ]
                                    }
                        }

def test_gbk_no_gff_test():
    """Test genomic mode without gff as input.
    """
    print("*** Test genomic mode without gff as input ***")
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


def test_gbk_from_gff_test():
    """Test genomic mode with a gff as input.
    """
    print("*** Test genomic mode with gff as input ***")
    gbk_test = 'test_gff.gbk'

    gbk_creation(genome=FNA_INPUT,
                proteome=FAA_INPUT,
                annot=ANNOT_INPUT,
                gff=GFF_INPUT,
                org=ORG_NAME,
                gbk=gbk_test,
                gobasic='go-basic.obo')

    compare_two_gbks(EXPECTED_GBK_NO_GFF, gbk_test)
    os.remove(gbk_test)
    return

def test_gbk_from_gff_metagenomic_mode():
    """Test metagenomic mode.
    """
    print("*** Test metagenomic mode ***")
    gbk_dir_test = 'gbk_mg/'
    os.makedirs(gbk_dir_test)
    gbk_creation(genome=FNA_DIR,
                proteome=FAA_DIR,
                annot=ANNOT_INPUT,
                org=ORG_FILE,
                gff=GFF_DIR,
                gbk=gbk_dir_test,
                gobasic='go-basic.obo',
                dirmode=True,
                metagenomic_mode=True)

    check_gbks_from_dir(gbk_dir_test)

    shutil.rmtree(gbk_dir_test)
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

def test_gbk_from_dir_test():
    """Test genomic mode with directories.
    """
    print("*** Test genomic mode with input directories without gff ***")
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

    check_gbks_from_dir(gbk_dir_test)

    shutil.rmtree(gbk_dir_test)
    return

def test_gbk_metagenomic_mode():
    """Test metagenomic mode.
    """
    print("*** Test metagenomic mode ***")
    gbk_dir_test = 'gbk_mg/'
    os.makedirs(gbk_dir_test)
    gbk_creation(genome=FNA_DIR,
                proteome=FAA_DIR,
                annot=ANNOT_INPUT,
                org=ORG_FILE,
                gff=None,
                gbk=gbk_dir_test,
                gobasic='go-basic.obo',
                dirmode=True,
                metagenomic_mode=True)

    check_gbks_from_dir(gbk_dir_test)

    shutil.rmtree(gbk_dir_test)
    return

def check_gbks_from_dir(gbk_dir):
    """Check if annotations in each gbk file are consistent with the expected ones.

    Args:
        gbk_dir (str): path to gbk directory
    """
    for gbk in os.listdir(gbk_dir):
        loaded_gbk = SeqIO.to_dict(SeqIO.parse(f"{gbk_dir}/{gbk}", "genbank"))
        annotations = {i:None 
                        for i in loaded_gbk.keys()}
        for gene in annotations: 
            for index in range(0, len(loaded_gbk[gene].features)):
                if loaded_gbk[gene].features[index].type=='CDS':
                    annotations[gene] = loaded_gbk[gene].features[index].qualifiers
                    # check annotations
                    for ann in ANNOTATIONS_TYPES:
                        if ann in annotations[gene]:
                            # print(set(annotations[gene][ann]), set(ANNOTATIONS_BY_GENOME[gene][ann]))
                            assert set(annotations[gene][ann]) == set(ANNOTATIONS_BY_GENOME[gene][ann])
    return

    
if __name__ == "__main__":
    test_gbk_no_gff_test()
    # test_gbk_from_gff_test()
    #test_gbk_from_dir()
    #test_gbk_metagenomic_mode()
