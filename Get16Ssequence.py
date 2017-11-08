

"""
##gff-version 3
NODE_10_length_3295_cov_425.534	barrnap:0.7	rRNA	5	114	1.8e-14	-	.	Name=5S_rRNA;product=5S ribosomal RNA
NODE_10_length_3295_cov_425.534	barrnap:0.7	rRNA	197	3124	0	-	.	Name=23S_rRNA;product=23S ribosomal RNA
NODE_12_length_1324_cov_73.664	barrnap:0.7	rRNA	1252	1323	6.5e-07	-	.	Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 60 percent of the 5S ribosomal RNA
NODE_16_length_898_cov_456.518	barrnap:0.7	rRNA	1	897	1.4e-262	-	.	Name=16S_rRNA;product=16S ribosomal RNA (partial);note=aligned only 56 percent of the 16S ribosomal RNA
NODE_22_length_586_cov_56.7485	barrnap:0.7	rRNA	514	585	6.5e-07	-	.	Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 60 percent of the 5S ribosomal RNA
NODE_4_length_267895_cov_61.0199	barrnap:0.7	rRNA	267823	267894	6.5e-07	-	.	Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 60 percent of the 5S ribosomal RNA
NODE_5_length_154291_cov_61.9736	barrnap:0.7	rRNA	2	73	6.5e-07	+	.	Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 60 percent of the 5S ribosomal RNA
NODE_7_length_109377_cov_57.9658	barrnap:0.7	rRNA	2	73	6.5e-07	+	.	Name=5S_rRNA;product=5S ribosomal RNA (partial);note=aligned only 60 percent of the 5S ribosomal RNA


"""
import os
import sys
import string
import argparse

def revcomplement(sequence):
	""" Reverse Complement """
	basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
	nucleotides = list(sequence)
	nucleotides = [basecomplement[base] for base in nucleotides]
	return ''.join(nucleotides)
	
	return revcomp

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

### END OF CLASS

class RibosomalGetter:
	
	def __init__(self, scaffoldsfile, barnappfile):
		self.scaffoldsfile = scaffoldsfile
		self.outfile = scaffoldsfile.rsplit('.', 1)[0] + '.16S.fasta'
		self.barnappfile = barnappfile
		
	def run(self):
		""" """
		riboDict = {}
		contiglist = []
		with open(self.barnappfile) as gfffile:
			for line in gfffile.readlines():
				if "product=16S" in line:
					columns = line.strip().split('\t')
					contigid = columns[0]
					start = int(columns[3])
					stop = int(columns[4])
					strand = columns[6]
					riboDict.setdefault(contigid, []).append((start, stop, strand))
					#print riboDict
		riboOut = self._acquireContigs(riboDict)
		FastaWriter(riboOut, 100, self.outfile).writedown() 
		return
		
		
		
	def _acquireContigs(self, riboDict):
		""" """
		riboOut = {}
		FS = FastaStream(self.scaffoldsfile)
		for fasta in FS:
			title = fasta.split('\n', 1)[0].strip().replace('>', '')
			if title in riboDict:
				sequence = "".join(fasta.split('\n')[1:])
				for coord in riboDict[title]:
					riboseq = sequence[coord[0]-1: coord[1]]
					if coord[2] == "-":
						riboseq = revcomplement(riboseq.upper())
					riboOut[title+'_'+str(coord[0])+'_'+str(coord[1])+'_'+coord[2]] = riboseq
		return riboOut
		

### END OF CLASS

def main(argv):
    """ main function """
    parser = argparse.ArgumentParser()

    parser.add_argument('-s', dest='scaffoldsfile', action='store', required=True, help='Consensus File Mandatory ')
    parser.add_argument('-b', dest='barnappfile', action='store', required=True, help='expected GenomeSize ')
    
    args = parser.parse_args()

    S = RibosomalGetter(args.scaffoldsfile, args.barnappfile)
    S.run()
    return


################
if __name__ == "__main__":
    sys.exit(main(sys.argv))
