#!/usr/bin/env python2.7



from sys import argv
import csv
from collections import defaultdict
import networkx
from networkx.algorithms.components.connected import connected_components

blast_in = argv[1]
query_in = argv[2]

query_lengths = {}
with open(query_in) as fin:
    for line in fin.readlines():
        line = [l.strip() for l in line.split("\t")]
        query_lengths[line[0]] = int(line[1])


outfmt6 = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


blast_table_filtered = []
connections = []
clustered = []
singletons = set([])
with open(blast_in) as fin, open(blast_in+".filtered", "w") as fout:
    blast_table = csv.DictReader(fin,fieldnames=outfmt6,dialect='excel-tab')
    blast_table_filtered = csv.DictWriter(fout, fieldnames=outfmt6,dialect='excel-tab')
    blast_table_filtered.writeheader()
    for row in blast_table:
        if row["qseqid"] != row["sseqid"]:

            query_length = query_lengths[row["qseqid"]]

            if int(row["length"]) >= 0.8*query_length and float(row["pident"]) >= 80:
                blast_table_filtered.writerow(row)
                connections.append((row["qseqid"],row["sseqid"]))
                clustered.append(row["qseqid"])
                clustered.append(row["sseqid"])

# second itertion
with open(blast_in) as fin, open(blast_in + ".filtered", "w") as fout:
    blast_table = csv.DictReader(fin, fieldnames=outfmt6, dialect='excel-tab')
    for row in blast_table:
        if row["qseqid"] == row["sseqid"]:
            if row["qseqid"] not in clustered:
                singletons.add(row["qseqid"])

with open(blast_in+".singletons", "w") as singletons_fin:
    singletons_fin.write("\n".join(singletons))

# create Graph of pairs
def to_graph(CL, edge):
    G = networkx.Graph()
    for part in CL:
        G.add_nodes_from(part)
        G.add_edges_from(edge)
    return G


G = to_graph(connections, connections)
connection_list = [list(se) for se in connected_components(G)]
connection_list.sort(key=len, reverse=True)
for se_list in connection_list:
    print ".fa ".join(se_list) + ".fa "


