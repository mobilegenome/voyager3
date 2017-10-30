#!/usr/bin/env python2


class SeqAnalyses:

    def __init__(self, seqid, fasta_in):
        from Bio import SeqIO

        from os import makedirs

        self.fasta_in = fasta_in
        self.seq = SeqIO.read(fasta_in, "fasta")
        self.outdir = seqid
        try:
            makedirs(seqid)
        except:
            pass

    @property
    def repmask(self):
        import subprocess
        cmdl = "RepeatMasker -species mammal %s" %self.fasta_in
        print(" ".join(cmdl.split(" ")))
        subprocess.check_call(cmdl.split(" "))

        return 1

    def blastNR(self):
        import httplib
        httplib.HTTPConnection._http_vsn = 10
        httplib.HTTPConnection._http_vsn_str = 'HTTP/1.0'
        from Bio.Blast import NCBIWWW
        print "Blasting"
        print self.seq.format("fasta")
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



#class CompAnalyses:
#