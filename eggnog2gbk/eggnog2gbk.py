"""Main module."""

from eggnog2gbk.gff_to_gbk import gff_to_gbk
from eggnog2gbk.fa_to_gbk import faa_to_gbk

def gbk_creation(genome:str, proteome:str, annot:str, gff:str, org:str, gbk:str, gobasic:str):
    if gff:
        gff_to_gbk(genome, proteome, annot, gff, org, gbk, gobasic)
    else:
        faa_to_gbk(genome, proteome, annot, org, gbk, gobasic)