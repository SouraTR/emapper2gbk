"""Main module."""

from eggnog2gbk import gff_to_gbk
from eggnog2gbk import fa_to_gbk

def gbk_creation(genome:str, proteome:str, annot:str, gff:str, org:str, gbk:str, gobasic:str):
    if gff:
        gff_to_gbk.main(genome, proteome, annot, gff, org, gbk, gobasic)
    else:
        fa_to_gbk.main(genome, proteome, annot, org, gbk, gobasic)