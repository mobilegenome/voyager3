#!/usr/bin/env python2

import os
import re
import signal
from Bio import SeqIO

print "load RepeatMasker library"
rmlib = SeqIO.index("/home/fritjof/programs/RepeatMasker/Libraries/RepeatMasker.lib", "fasta")
print "done."


class SeqAnalyses:
    def __init__(self, seqid, fasta_in):
        from Bio import SeqIO

        from os import makedirs

        self.fasta_in = fasta_in
        self.seqid = seqid
        self.RMhitfiles = []

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
                    rmhitname = hit.group(1)
                    rmhitsummary = "{elem}:{start}-{end}({score})".format(elem=rmhitname,
                                                                 start=feature["start"],
                                                                 end=feature["end"],
                                                                 score=100 - feature["score"])
                    rmhits.append(rmhitsummary)

                    rmlib_hits = (k for k in rmlib.keys() if re.match(rmhitname, k))
                    for hit in rmlib_hits:
                        SeqIO.write(rmlib[hit], "%s/RM.%s.fa" % (self.seqid, rmhitname), "fasta")
                        self.RMhitfiles.append("%s/RM.%s.fa" % (self.seqid, rmhitname))

                    rmhits.append(rmhitsummary)
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
        print("Blasting")
        try:
  #          signal.signal(signal.SIGALRM, handler)
   #         signal.alarm(60)
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
            print("BLAST error. Continuing")
            return "NA"
        # except "IncompleteReadError":
        #     print("BLAST timeout. Continuing")
        #     return "NA"
        # except TimeoutError:
        #     print("BLAST timeout. Continuing")
        #     return "NA"


    def dotplot(self):
        window = 11
        seqfile1 = self.fasta_in
        seq1 = SeqIO.read(seqfile1, "fasta")
        seq2 = SeqIO.read(self.RMhitfiles[0], "fasta")

        seq_one = str(seq1.seq).upper()
        seq_two = str(seq2.seq).upper()
        data = [[(seq_one[i:i + window] <> seq_two[j:j + window]) for j in range(len(seq_one) - window)] for i in range(len(seq_two) - window)]

        import pylab
        pylab.gray()
        pylab.imshow(data)
        pylab.xlabel("%s (length %i bp)" % (seq1.id, len(seq1)))
        pylab.ylabel("%s (length %i bp)" % (seq2.id, len(seq2)))
        pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
        pylab.show()

    def call_dotplot(self):

        seqfile1 = self.fasta_in
        seq1 = SeqIO.read(seqfile1, "fasta")
        if len(set(self.RMhitfiles)) ==1:
            seq2 = SeqIO.read(self.RMhitfiles[0], "fasta")
            self.dotplot2(seq1, seq2)
        elif 1 < len(set(self.RMhitfiles)) < 5:
            for rmseqfile in set(self.RMhitfiles):
                seq2 = SeqIO.read(rmseqfile, "fasta")
                self.dotplot2(seq1, seq2)
        else:
            return 0

    def dotplot2(self, seq1, seq2):

        window = 7
        dict_one = {}
        dict_two = {}

        for (seq, section_dict) in [(str(seq1.seq).upper(), dict_one),
                                    (str(seq2.seq).upper(), dict_two)]:
            for i in range(len(seq) - window):
                section = seq[i:i + window]
                try:
                    section_dict[section].append(i)
                except KeyError:
                    section_dict[section] = [i]

        dict_one_rev = {}
        dict_two_rev = {}

        for (seq, section_dict) in [(str(seq1.seq).upper(), dict_one_rev),
                                    (str(seq2.seq.reverse_complement()).upper(), dict_two_rev)]:
            for i in range(len(seq) - window):
                section = seq[i:i + window]
                try:
                    section_dict[section].append(i)
                except KeyError:
                    section_dict[section] = [i]

        matches = set(dict_one).intersection(dict_two)
        matches_rev = set(dict_one_rev).intersection(dict_two_rev)
        print("%i unique matches" % len(matches))

            # Create lists of x and y co-ordinates for scatter plot
        x = []
        y = []
        for section in matches:
            for i in dict_one[section]:
                for j in dict_two[section]:
                    x.append(i)
                    y.append(j)

        xr = []
        yr = []
        for section in matches_rev:
            for i in dict_one_rev[section]:
                for j in dict_two_rev[section]:
                    xr.append(i)
                    yr.append(j)

        import pylab
        pylab.cla()  # clear any prior graph
        pylab.gray()
        pylab.scatter(x, y, c="green", s=0.5)
        pylab.scatter(xr, yr, c="red", s=0.5)
        pylab.xlim(0, len(seq1) - window)
        pylab.ylim(0, len(seq2) - window)
        pylab.xlabel("%s (length %i bp)" % (seq1.id, len(seq1)))
        pylab.ylabel("%s (length %i bp)" % (seq2.id, len(seq2)))
        pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
        pylab.savefig(os.path.join(self.outdir, "%s.%s.dotplot.png" %(seq1.id, seq2.id.split("#")[0])))
        return 1

    def secstruct(self):

        return 1

def handler(signum, frame):
    print "Timeout."
    raise Exception("TimeoutError")

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
