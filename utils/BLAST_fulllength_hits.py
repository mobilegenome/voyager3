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



def to_graph(CL, edge):
    G = networkx.Graph()
    for part in CL:
        G.add_nodes_from(part)
        G.add_edges_from(edge)
    return G


G = to_graph(connections, connections)
for se in connected_components(G):
    print list(se)



# conn_dict = defaultdict(list)
# for tup in connections:
#     q = tup[0]
#     s = tup[1]
#
#     if q in conn_dict.keys():
#         conn_dict[q].append(s)
#     elif s in conn_dict.keys():
#         conn_dict[s].append(q)
#     else:
#         if q in conn_dict.values():
#             q_hits = [k for k,v in conn_dict.iteritems() if q in v]
#             assert (len(q_hits) == 1)
#             conn_dict[q_hits[0]].append(s)
#         elif s in conn_dict.values():
#             s_hits = [k for k, v in conn_dict.iteritems() if s in v]
#             assert (len(q_hits) == 1)
#             conn_dict[s_hits[0]].append(q)
#         else:
#             conn_dict[q].append(s)
# q
#
#
#
# for k,v in conn_dict.iteritems():
#     print "%s\t%s" %(k, ",".join(v))
# q
