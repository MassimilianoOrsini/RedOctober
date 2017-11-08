"""

"""

import os
import sys
import string
import math
import time
import re
import argparse

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





class CheckContigs:
	
	def __init__(self, infilename,genomesize, mode, label=""):
		self.infilename = infilename.strip()
		self.genomesize = genomesize
		self.mode = mode
		self.label = label
		self.contigs = []
		self.assembled = []
		self.headline = ["# Sample", "Len.", "Assembled", "AssembledPerc.", "Contigs", "Avg.", "Median", "Min.", "Max.", "Lt200", "N50", "Gaps", "GapsAvg", "GapsMed", "GapsMax", "GapsSt20"]
		
	def analyze(self):
		""" """
		self._acquire()
		self._countsLinear()
		return

		
	def _acquire(self):
		""" """
		fastaStream = FastaStream(self.infilename)
		self.checkcounter = 0
		for fasta in fastaStream:
			head = fasta.split("\n")[0]
			seq =  "".join(fasta.split("\n")[1:])
			self.contigs.append(len(seq))
			if self.mode == "atgc":
				self.assembled.append((seq.upper().count("A")+seq.upper().count("C")+seq.upper().count("T")+seq.upper().count("G")))
			else:
				self.assembled.append(len(seq) - seq.upper().count("N"))
			self.checkcounter += 1
		return

	def _countsLinear(self):
		""" """
		print "\t".join(self.headline)
		if self.checkcounter != 0:
			resline = ["# "+self.infilename.split("/")[-1]+self.label, str(self._basesAssembled(self.contigs)), str(self._basesAssembled(self.assembled)), \
			str(round((float(self._basesAssembled(self.assembled))/float(self._basesAssembled(self.contigs))*100), 1)), str(len(self.contigs)), \
			str(self._calcAverage(self.contigs)) , str(self._calcMedian(self.contigs)), str(min(self.contigs)), str(max(self.contigs)), \
			str(self._calcOver200(self.contigs)), str(self._calcN50(self.contigs)), "NA", "NA", "NA", "NA", "NA"]
			
		else:
			resline = ["# "+self.infilename.split("/")[-1], "0", "0", "0", "0", "0" , "0", "0", "0", "0", "0", "NA", "NA", "NA", "NA", "NA"]
		print "\t".join(resline)
		return

	def _calcAverage(self, values):
		""" """
		if len(values) == 0: return "-"
		total = 0
		for element in values:
			total += element
		average = total/float(len(values))
		return int(round((average), 1))
		
	def _calcMedian(self, values):
		""" """
		if len(values) == 0: return "-"
		theValues = sorted(values)
		if len(theValues) % 2 == 1:
			return int(theValues[(len(theValues)+1)/2-1])
		else:
			lower = theValues[len(theValues)/2-1]
			upper = theValues[len(theValues)/2]
			return int((float(lower + upper)) / 2)	

	def _calcN50(self, values):
		""" """
		values.sort()
		theValues = []
		for value in values:
			theValues += [value] * value
		if len(theValues) % 2 == 0:
			medianpos = len(theValues)/2
			N50 = int(float(theValues[medianpos] + theValues[medianpos-1]) /2)
		else:
			medianpos = len(theValues)/2
			N50 = theValues[medianpos]
		return N50

	def _calcNG50(self, values):
		""" """
		theValues = sorted(values)
		theValues.reverse()
			
		index = 0
		tmpG50 = 0
		try:
			while tmpG50 < float(self.genomesize)/2:
				tmpG50 += theValues[index]
				index += 1
			NG50 =  theValues[index]
		except:
			NG50 = 0
		return NG50

	def _basesAssembled(self, values):
		""" """
		total = 0
		for element in values:
			total += element
		return total

	def _calcOver200(self, values):
		""" """
		total = 0
		for element in values:
			if element >= 200: total +=1
		return total

	def _calcOver2000(self, values):
		""" """
		total = 0
		for element in values:
			if element >= 2000: total +=1
		return total

	def _calcGaps5(self, values):
		""" """
		total = 0
		for element in values:
			if element <= 5: total +=1
		return total

### END OF CLASS

def main(argv):
    """ main function """
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', dest='fasta', action='store', required=True, help='Consensus File Mandatory ')
    parser.add_argument('-g', dest='genomesize', action='store', required=False, default='3000000', help='expected GenomeSize ')
    parser.add_argument('-m', dest='mode', action='store', required=False, help='report file ')
    parser.add_argument('-l', dest='label', action='store', required=False, default='', help='sample label ')
    
    args = parser.parse_args()

    S = CheckContigs(args.fasta, args.genomesize, args.mode, args.label)
    S.analyze()
    return


################
if __name__ == "__main__":
    sys.exit(main(sys.argv))
