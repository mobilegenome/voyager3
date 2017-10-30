#!/usr/bin/env python2.7

import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from sys import argv

infile = argv[1]
taxon = argv[2]

with open(infile) as fin:
    in_csv = csv.DictReader(fin)

    for row in in_csv:
        seq = Seq(row["Consensus"])
        record = SeqRecord(seq,
                           id=row["Cluster"],
                           name=row["Cluster"],
                           description="length: {length} bp, genome_proportion: {prop:.5}%, GC: {gc:.4}% ".format(
                               length=float(row["Consensus length"]), prop=float(row["Genome Proportion[%]"]),
                               gc=GC(seq)))

        with open("%s.%s.fa" % (taxon, row["Cluster"]), "w") as fout:
            fout.write(">%s %s\n" % (record.id, record.description))
            fout.write(str(record.seq) + "\n")
