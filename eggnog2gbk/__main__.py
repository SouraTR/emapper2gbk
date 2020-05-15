"""Console script for eggnog2gbk."""
import argparse
import sys
import pkg_resources
from eggnog2gbk.eggnog2gbk import gbk_creation

VERSION = pkg_resources.get_distribution("eggnog2gbk").version
LICENSE = """Copyright (C) Pleiade and Dyliss Inria projects
MIT License - Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions. See LICENSE for more details \n
"""
MESSAGE = """
Starting from fasta and eggnog annotation files, build a gbk file that is suitable for metabolic network reconstruction with Pathway Tools. Adds the GO terms and EC numbers annotations in the genbank file.
"""


def cli():
    """Console script for eggnog2gbk."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + VERSION + "\n" + MESSAGE)

    parser.add_argument(
        "-fp",
        "--fastaprot",
        help="faa file",
        required=True,
        type=str)

    parser.add_argument(
        "-fg",
        "--fastagenome",
        help="fna file",
        required=True,
        type=str)

    parser.add_argument(
        "-g",
        "--gff",
        help="gff file",
        required=True,
        type=str)

    parser.add_argument(
        "-a",
        "--annotation",
        help="eggnog annotation file",
        required=True,
        type=str)

    parser.add_argument(
        "-o",
        "--output",
        help="gbk output file",
        required=True,
        type=str)

    parser.add_argument(
        "-go",
        "--gobasic",
        help="go ontology, will be downloaded if not provided",
        required=False,
        type=str)

    parser.add_argument(
        "-n",
        "--name",
        help="organism/genome name",
        required=True,
        type=str)

    args = parser.parse_args()

    gbk_creation(args.fastagenome, args.fastaprot, args.annotation, args.gff, args.name, args.output)

if __name__ == "__main__":
    sys.exit(cli())  
