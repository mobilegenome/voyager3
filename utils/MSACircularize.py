#!/usr/bin/env python

'''
MSACircularize.py
(c) Fritjof Lammers
Trims  multiple sequence alignments (MSA) of circular sequences like mtDNA, plasmid, or satDNA monomer.
The algorithms double each sequence and subsequently strips it down to its original lengths.
Inner-sequence gaps are maintained, but might require re-alignment of the processed MSA with external tools.

Handles FASTA files.

Usage:
python MSACircularize.py INPUTFILE OUTPUTFILE
'''

from sys import argv
import re

from Bio import AlignIO
from Bio import Seq

scriptname = argv[0].split("/")[-1]
infile = argv[1]
outfile = argv[2]

msa_format = "fasta"
aln_metadata = []

print "[ %s ] Open Alignment... " % scriptname,
with open(infile) as fin:
    print "process... "
    aln = AlignIO.read(fin, msa_format)  # load MSA

    #  iteration to save gap-states in alignment
    gs_max = 0
    for record in aln:

        gapmatch = re.match("[-]+", str(record.seq))
        gs = len(gapmatch.group(0)) if gapmatch else 0
        gs_max = gs if gs_max < gs else gs_max

    # second teration to alter sequences in alignment

    for record in aln:
        seqlen = len(str(record.seq).strip("-"))  # original sequence length
        record.seq = Seq.Seq(str(record.seq).rstrip("-") + str(record.seq).strip("-"))  # double sequence to circularize
        record.seq = record.seq[gs_max:(seqlen + gs_max)]


with open(outfile, "w") as fout:
    print "[ %s ] Save alignment to %s ..." % (scriptname, outfile)
    AlignIO.write(aln, fout, msa_format)
    print "[ %s ] Done." % scriptname
