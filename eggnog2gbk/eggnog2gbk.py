"""Main module."""

import gffutils
import numpy as np
import os
import pandas as pa
import pronto
import re
import requests
import shutil

from Bio import SeqFeature as sf
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

from eggnog2gbk.gff_to_gbk import gff_to_gbk

def gbk_creation(genome:str, proteome:str, annot:str, gff:str, org:str, gbk:str):
    gff_to_gbk(genome, proteome, annot, gff, org, gbk)