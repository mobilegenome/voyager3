#!/usr/bin/env python2

import os
import re
import signal


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
                    # TODO extract and save repeat sequence
                    rmhits.append(hit)
                return ",".join(rmhits)

        except "CalledProcessError":
            print("RepeatMasker exited  with error")
            return 0

    def blastNR(self):
        import httplib
        from Bio.Blast import NCBIWWW
        from Bio.Blast import NCBIXML
        httplib.HTTPConnection._http_vsn = 10
        httplib.HTTPConnection._http_vsn_str = 'HTTP/1.0'
        print
        "Blasting"
        print
        self.seq.format("fasta")
        try:
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(60)
            result_handle = NCBIWWW.qblast("blastn", "nt", self.seq.format("fasta"))
            blast_records = NCBIXML.parse(result_handle)
            if blast_records:
                blast_record = blast_records.next()
                alignment =  blast_record.alignments[0]
                for hsp in alignment.hsps:
                    if hsp.expect < 0.004:
                        print('****Alignment****')
                        print('sequence:', alignment.title)
                        print('length:', alignment.length)
                        print('e value:', hsp.expect)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')
                        firstblasthit = "{acc},{alnlength}".format(acc=alignment.title,
                                                                 alnlength=alignment.length)
                with open("%s/%s.blast.xml" % (self.seqid, os.path.basename(self.fasta_in)),
                          "w") as blastout_file:
                    blastout_file.write(result_handle.read())
                    print("BLASTn results written to %s/%s.blast.xml" %(self.seqid, os.path.basename(self.fasta_in)))
                    result_handle.close()
                return firstblasthit
            else:
                return "NA"
        except:
            print("BLAST failed. Continuing")
            return "NA"


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

def handler(signum, frame):
    print "Forever is over!"
    raise Exception("end of time")

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
