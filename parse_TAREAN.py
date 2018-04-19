#!/usr/bin/env python2

import csv
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from sys import argv
from y3 import SeqAnalyses


def csv_parser(infile, outfile):
    fin = open(infile)
    fout = open(outfile, "w")
    taxon = infile.split(".")[0]
    in_csv = csv.DictReader(fin)

    in_csv.fieldnames.append("GC%")
    in_csv.fieldnames.append("RM_hits")
    in_csv.fieldnames.append("BLASTn_hits(nr)")

    out_csv = csv.DictWriter(fout, fieldnames=in_csv.fieldnames)
    out_csv.writeheader()

    for row in in_csv:

        if len(row) == 0:
            continue
        if not row["Consensus"]:
            continue
        seq = Seq(row["Consensus"])
        row["GC%"] = "{gc:.4}".format(gc=GC(seq))
        try:
            record = SeqRecord(seq,
                               id=row["Cluster"],
                               name=row["Cluster"],
                               description="length: {length} bp, genome_proportion: {prop:.5}%, GC: {gc:.4}% ".format(
                                   length=float(row["Consensus length"]), prop=float(row["Genome Proportion[%]"]),
                                   gc=row["GC%"]))
        except KeyError:
            record = SeqRecord(seq,
                               id=row["Cluster"],
                               name=row["Cluster"],
                               description="length: {length} bp, genome_proportion: {prop:.5}%, GC: {gc:.4}% ".format(
                                   length=float(row["Consensuslength"]), prop=float(row["GenomeProportion[%]"]),
                                   gc=row["GC%"]))


        fasta_out = "%s.%s.fa" % (taxon, row["Cluster"])
        with open(fasta_out, "w") as fout:
            fout.write(">%s.%s %s\n" % (taxon, record.id, record.description))
            fout.write(str(record.seq)+"\n")

        seqinspect = SeqAnalyses("%s.%s" % (taxon, row["Cluster"]), fasta_out, species)
        row["RM_hits"] = seqinspect.repmask
        seqinspect.call_dotplot()
        row["BLASTn_hits(nr)"] = seqinspect.blastNR()

        out_csv.writerow(row)

        # TODO create HTML/PDF output
    fin.close()
    fout.close()

if __name__ == "__main__":
    infile = argv[1]
    species = argv[2]
    outfile = argv[3]
    csv_parser(infile, outfile)
