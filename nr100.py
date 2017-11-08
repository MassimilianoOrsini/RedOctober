import os
import string
import sys
import argparse
from operator import *

class FastaStream:

    def __init__(self, multifastafile):
        """ Generate a stream of 'FASTA strings' from an io stream. """
        self.infile = open(multifastafile)
        self.buffer = []

    def __del__(self) :
        self.infile.close()

    def __iter__(self):
        return self

    def next(self):
        while 1:
            try:
                line = self.infile.next()
                if line.startswith('>') and self.buffer:
                    fasta = "".join(self.buffer)
                    self.buffer = [line]
                    return fasta
                else:
                    self.buffer.append(line)
            except StopIteration:
                if self.buffer:
                    fasta = "".join(self.buffer)
                    self.buffer = []
                    return fasta
                else:
                    raise StopIteration

### END OF CLASS

class FastaWriter:

    def __init__(self, fastasdict, linewidth, outname, mode="w"):
        """ writes fasta dicts in files """
        self.fastasdict = fastasdict
        self.linewidth = int(linewidth)
        self.outname = outname
        self.mode = mode

    def writedown(self):
        """ main method """
        with open(self.outname, self.mode) as outfile:
            for fastatitle in self.fastasdict:
                sequence = self.fastasdict[fastatitle]
                if not fastatitle.startswith(">"):
                    fastatitle = ">"+fastatitle
                outfile.write(fastatitle+ '\n')
                tmpline = ""
                for nucleotide in sequence:
                    tmpline += nucleotide
                    if len(tmpline) == self.linewidth:
                        outfile.write(tmpline + '\n')
                        tmpline = ""
                if tmpline != "": outfile.write(tmpline + '\n')
        return

### END OF CLASS



class NonRedundantFasta:
    
    def __init__(self, infile):
        self.infile = infile
        self.outfile = infile.rsplit('.', 1)[0] +'.nr.fas'
        
    def nr(self):
        """ """
        nrdict = {}
        FS = FastaStream(self.infile)
        for fasta in FS:
            title = fasta.split('\n')[0]
            sequence = "".join(fasta.split('\n')[1:]).upper()
            if sequence not in nrdict:
                contained = False
                for element in nrdict:
                    if sequence in element:
                        contained = True
                if contained == False:
                    nrdict[sequence] = title
        trasposeddict = {}
        counter = 1
        for sequence, title in nrdict.items():
            #trasposeddict[title] = sequence
            trasposeddict["C"+str(counter)] = sequence
            counter += 1
        FastaWriter(trasposeddict, 100, self.outfile).writedown()
        return
    
    
    
if __name__ == "__main__":
    infile = sys.argv[1]
    NR = NonRedundantFasta(infile)
    NR.nr()