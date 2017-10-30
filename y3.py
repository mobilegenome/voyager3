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

    def repmask(self):
        import subprocess
        cmdl = "RepeatMasker -species mammal %s " %self.fasta_in
        subprocess.run(cmdl, shell=True, check=True)

        return 1

    def blastNR(self):
        from Bio.Blast import NCBIWWW
        print "Blasting"
        print self.seq.format("fasta")
        result_handle = NCBIWWW.qblast("blastn", "nr", self.seq.format("fasta"))

        with open("my_blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
            result_handle.close()
        return 1

    def dotplot(self):

        return 1

    def secstruct(self):

        return 1



#class CompAnalyses:
#