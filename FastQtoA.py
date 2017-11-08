import os
import string
import sys
import argparse



class FastqStream:

	def __init__(self, multifastqfile):
		"""
		Generate a stream of 'FASTQ strings' from an io stream.
		"""
		self.infile = open(multifastqfile)
		self.buffer = []

	def __del__(self) :
		self.infile.close()

	def __iter__(self):
		return self

	def next(self):
		while 1:
			try:
				head = self.infile.next()
				seq = self.infile.next()
				comment = self.infile.next()
				qual = self.infile.next()
				self.buffer = [head, seq, comment, qual]
				fastq = "".join(self.buffer)
				return fastq
			except StopIteration:
				if self.buffer:
					fastq = "".join(self.buffer)
					self.buffer = []
					return fastq
				else:
					raise StopIteration


############################


class Converter:
	
	def __init__(self, infilename, outname=None):
		self.infilename = infilename
		if not outname:
			self.outname = infilename.rsplit('.', 1)[0] + '.fa'
		else:
			self.outname = outname
		
	def convert(self):
		""" """
		multifastq = FastqStream(self.infilename)
		with open(self.outname, 'wb') as out:
			for fastq in multifastq:
				title = ">"+fastq.split('\n')[0].replace("@", "")
				sequence = fastq.split('\n')[1].strip()
				out.write(title + '\n' + sequence +'\n')
		return


# END OF CLASS



def main(argv):
	# main function
	parser = argparse.ArgumentParser()

	parser.add_argument('-i', dest='infilename', action='store', required=True, help='Draft File Mandatory ')
	parser.add_argument('-o', dest='outname', action='store', required=False, default=None, help='Vcf File Mandatory ')
	
	args = parser.parse_args()

	FG = Converter(args.infilename, args.outname)
	FG.convert()
	return

################
if __name__ == "__main__":
	sys.exit(main(sys.argv))	 
	 
