#!/usr/bin/env python2

import os
import re


class SeqAnalyses:
    def __init__(self, seqid, fasta_in):
        from Bio import SeqIO

        from os import makedirs

        self.fasta_in = fasta_in
        self.seqid = seqid

        self.seq = SeqIO.read(fasta_in, "fasta")
        self.outdir = seqid
        try:
            makedirs(seqid)
        except:
            pass

    @property
    def repmask(self):
        import subprocess
        cmdl = "RepeatMasker -species mammal -dir {outdir} -gff -html {seqfa}".format(outdir=self.seqid,
                                                                                      seqfa=self.fasta_in)
        print("Searching for repeats in %s" % self.seqid),
        print(" ".join(cmdl.split(" ")))

        try:
            rm_stdout = subprocess.check_output(cmdl.split(" "))
            if rm_stdout.splitlines()[-1].rstrip().startswith("No repetitive sequences were detected in"):
                print("RepeatMasker did not identify any repeats.")
            else:
                rmhits = []
                rmout_gff = "%s/%s.out.gff" % (self.seqid, os.path.basename(self.fasta_in))
                gff = parseGFF(rmout_gff)
                for feature in gff:
                    elem_re = re.compile("\:(.*)\"")
                    hit = elem_re.search(feature["other"])
                    hit = "{elem}:{start}-{end}({score})".format(elem=hit.group(1),
                                                                 start=feature["start"],
                                                                 end=feature["end"],
                                                                 score=100 - feature["score"])
                    rmhits.append(hit)
                return ",".join(rmhits)

        except "CalledProcessError":
            print("RepeatMasker exited  with error")
            return 0

    def blastNR(self):
        import httplib
        httplib.HTTPConnection._http_vsn = 10
        httplib.HTTPConnection._http_vsn_str = 'HTTP/1.0'
        from Bio.Blast import NCBIWWW
        print
        "Blasting"
        print
        self.seq.format("fasta")
        try:
            result_handle = NCBIWWW.qblast("blastn", "nr", self.seq.format("fasta"))
        except IncompleteReadError:
            print("BLAST failed. Continuing")

        with open("my_blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
            result_handle.close()
        return 1

    def dotplot(self):
        # Create lists of x and y co-ordinates for scatter plot
        x = []
        y = []
        for section in matches:
            for i in dict_one[section]:
                for j in dict_two[section]:
                    x.append(i)
                    y.append(j)

        import pylab
        pylab.cla()  # clear any prior graph
        pylab.gray()
        pylab.scatter(x, y)
        pylab.xlim(0, len(rec_one) - window)
        pylab.ylim(0, len(rec_two) - window)
        pylab.xlabel("%s (length %i bp)" % (rec_one.id, len(rec_one)))
        pylab.ylabel("%s (length %i bp)" % (rec_two.id, len(rec_two)))
        pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
        pylab.show()

        return 1

    def secstruct(self):

        return 1


def parseGFF(infile):
    gff = []
    with open(infile) as fin:
        for line in fin.readlines():
            if line.startswith("#"):
                continue
            else:
                line = line.split("\t")
                feature = {"seqid": line[0],
                           "source": line[1],
                           "type": line[2],
                           "start": int(line[3]),
                           "end": int(line[4]),
                           "score": float(line[5]),
                           "strand": line[6],
                           "phase": line[7],
                           "other": line[8]}
                gff.append(feature)
    return gff


    # class CompAnalyses:
    #
