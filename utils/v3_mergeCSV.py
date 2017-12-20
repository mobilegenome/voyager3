#!/usr/bin/env python2.7

import csv
from sys import argv
from Bio.SeqUtils import GC
from Bio.Seq import Seq

csv_list = argv[1:]

def UnicodeDictReader(utf8_data, **kwargs):
    csv_reader = csv.DictReader(utf8_data, **kwargs)
    for row in csv_reader:
        yield {unicode(key, 'utf-8'): unicode(value, 'utf-8') for key, value in row.iteritems()}


# 21 -num,Cluster,GenomeProportion[%],GenomeProportionAdjusted[%],Sizereal,Satelliteprobability,Consensuslength,Consensus,
# Graphlayout,Kmeranalysis,Connectedcomponentindex C,PaircompletenessindexP,Kmercoverage,|V|,|E|,
# Pbsscore,Similaritybasedannotation,dataset,GC%,RM_hits,BLASTn_hits(nr)
# 21 - num,Cluster,Genome Proportion[%],Size real,Satellite probability,Consensus length,Consensus,
# Kmer analysis,Graph layout,Connected component index C,Pair completeness index P,Kmer coverage,|V|,|E|,
# Pbs score,The longest ORF length,Similarity based annotation,dataset,GC%,RM_hits,BLASTn_hits(nr)



common_fieldnames = "sample", "num", "Cluster", "GenomeProportion[%]", "Sizereal", "Satelliteprobability", \
                    "Consensuslength", "Consensus", "Similaritybasedannotation", "dataset", "GC%", "RM_hits", "BLASTn_hits(nr)"

out_csv_file = "all_out.tsv"
collected_rows = []
for in_csv_file in csv_list:
    with open(in_csv_file) as fin:
        sample = in_csv_file.split(".")[0]
        in_csv = csv.DictReader(fin)
        # remove spaces from fieldnames to make comparable
        in_csv.fieldnames = [elem.replace(" ", "") for elem in in_csv.fieldnames]

        for row in in_csv:
            row_selectedfields = {}
            for fieldname, cell in row.iteritems():
                # fieldname = unicode(fieldname, "utf-8")
                fieldname = fieldname.replace(" ", "")
                if fieldname in common_fieldnames:
                    if fieldname.strip() == "GC%" and not cell.strip():
                        row_selectedfields["GC%"] = "{gc:.4}".format(gc=GC(Seq(row["Consensus"])))
                        print
                        "GC% calculated and added"
                    elif cell.strip():
                        row_selectedfields[fieldname] =  cell
                    else:
                        row_selectedfields[fieldname] =  "NA"

            row_selectedfields.update({"sample":sample})
            collected_rows.append(row_selectedfields)


with open(out_csv_file, "w") as fout:
    out_csv = csv.DictWriter(fout, dialect="excel-tab", fieldnames=common_fieldnames)
    out_csv.writeheader()
    out_csv.writerows(collected_rows)


print "All Done."