import sys
import string
import os
import argparse


class FastQReaderMultiLine:

	def __init__(self, fastqfile):
		self.fastqfile = open(fastqfile)
				
	def next(self):
		""" """
		head = self.fastqfile.next().strip()
		sequence = ""
		seqlinecounter = 0
		row = self.fastqfile.next()
		while not row.startswith("+"):
			sequence += row
			seqlinecounter += 1
			row = self.fastqfile.next()
		comment = row #self.fastqfile.next()
		qual = ""
		qualcounter = 0
		while qualcounter != seqlinecounter:
			qual += self.fastqfile.next()
			qualcounter += 1
		return head, sequence, comment, qual

### END OF CLASS
		
class FastQtoAconverter:

	def __init__(self, fastqfile):
		self.fastqfile = fastqfile
		self.fastafile = fastqfile.rsplit(".", 1)[0] + ".fasta"
		
	def convert(self):
		""" """
		FQ = FastQReaderMultiLine(self.fastqfile)
		out = open(self.fastafile, "w")
		try:
			while True:
				fastq = FQ.next()
				head = fastq[0].strip().replace("@", ">").strip()
				sequence =  fastq[1].replace("\n", "").strip()
				out.write(head +"\n")
				tmpline = []
				for nucleotide in sequence:
					tmpline.append(nucleotide)
					if len(tmpline) >= 100:
						out.write("".join(tmpline) + "\n")
						tmpline = []
				out.write("".join(tmpline) + "\n")
		except StopIteration:
			out.close()
			return		
		
		
## END OF CLASS		
		
def main(argv):
    """ main function """
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', dest='fastqfile',  action='store', required=True,  help='FastQ file mandatory')
    
    args = parser.parse_args()
    FQ = FastQtoAconverter(args.fastqfile)
    FQ.convert()
    return


if __name__ == "__main__":
    sys.exit(main(sys.argv))	
    			
