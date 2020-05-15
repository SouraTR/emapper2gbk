#!/usr/bin/env python

"""Tests for `eggnog2gbk` package."""

import pytest
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

EXPECTED_GBK = "betaox_from_gff.gbk"
GBK_OUTPUT = "test.gbk"

def gbk_from_gff_test():
    """Test with a gff as input.
    """
    gbk_creation(genome=FNA_INPUT,
                proteome=FAA_INPUT,
                annot=ANNOT_INPUT,
                gff=GFF_INPUT,
                org=ORG_NAME,
                gbk=GBK_OUTPUT)

    loaded_gbk = SeqIO.parse(EXPECTED_GBK, "genbank")
    loaded_test = SeqIO.parse(GBK_OUTPUT, "genbank")

    #TODO finish

    return
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def gbk_no_gff_test():
    """Test without a gff as input.
    """
    #TODO
    return
