"""


"""


import os
import string
import sys
import time

class FastaStream:

    def __init__(self, multifastafile):
        """
        Generate a stream of 'FASTA strings' from an io stream.
        """
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

class FastAFilter:

    def __init__(self, infilename, keyword):
        self.infilename = infilename
        self.keyword = int(keyword)

    def filter(self):
        """ """
        with open(self.infilename+".filt.fasta", "w") as out:
            fastafile = FastaStream(self.infilename)
            for fasta in fastafile:
                head = fasta.split("\n")[0].strip()
                sequence = "".join(fasta.split("\n")[1:])
                if len(sequence) >= self.keyword:
                    out.write(fasta.strip()+'\n')
                






if __name__ == "__main__":
    infilename  = sys.argv[1]
    keyword = sys.argv[2]
    FTQ = FastAFilter(infilename, keyword)
    FTQ.filter()
    
