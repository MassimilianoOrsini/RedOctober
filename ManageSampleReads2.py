# -*- coding: utf-8 -*-
"""
This tool trims paired-end FASTQ files on the basis of quality score or left/right position, retaining mate integrity.
Reads without mate after filtering are saved in a separate output file.
"""

import math
import sys
import os
import string

def ht30(values):
	""" Arithmetic mean of a list of values """
	tmp = 0
	if values:
		for element in values:
			if element >= 30:
				tmp += 1
		return round((float(tmp) /float(len(values)) * 100), 1)
	else:
		return 0


def average(values):
	""" Arithmetic mean of a list of values """
	return float(sum(values)) / float(len(values)) if len(values) else float('nan')

def median(values):
	""" calc median over values """
	theValues = sorted(values)
	if len(theValues) % 2 == 1:
		return int(theValues[(len(theValues)+1)/2-1])
	else:
		try:
			lower = theValues[len(theValues)/2-1]
			upper = theValues[len(theValues)/2]
			return int((float(lower + upper)) / 2)
		except:
			lower = 0
			upper = 0
			return 0

def phred2sanger(phred_scores):
	""" Convert an array of Phred quality scores (integers in [0, 93]) to a Sanger-encoded quality string"""
	return ''.join([chr(score + 33) for score in phred_scores])


def sanger2phred(sanger_string):
	""" Convert a Sanger-encoded quality string to an array of Phred quality scores"""
	return [ord(ch) - 33 for ch in sanger_string]

############################


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

# END OF CLASS

class ManageSampleReads:
	
	
	def __init__(self, samplename):
		self.basename = samplename#.split("_resline.csv")[0]
		self.raw1 = samplename + "_R1_001.fastq" # in standby non lavora con i symlink
		self.raw2 = samplename + "_R2_001.fastq"
		self.trimm1 = samplename + "_R1.treads.fastq"
		self.trimm2 = samplename + "_R2.treads.fastq"
		self.unpaired = samplename + "_unp.treads.fastq"
		self.foreward = samplename + "_m.F.fastq"
		self.reverse = samplename + "_m.R.fastq"
		self.merged = samplename  + "_m.merg.fastq"
		self.working = samplename + "_working.fastq"
		self.passing = False
		self.readsfiles = [self.raw1, self.raw2, self.trimm1, self.trimm2, self.unpaired, self.foreward, self.reverse, self.merged, self.working] # in standby per i symlink
		self.checked = []
		
	def manage(self):
		""" """
		self._checkExistence()
		resHead = ['#Readsfile', 'NumOfReads', 'Assembled', 'Q30', 'AvgQual', 'MinLen', 'AvgLen', 'MedLen', 'MaxLen']
		with open(self.basename+'_samplereadscheck.tab', 'w') as outfile:
			outfile.write("\t".join(resHead) + "\n")
			print "\t".join(resHead)
			for readsfile in self.readsfiles:
				if readsfile in self.checked:
					numOfReads, assembled, q30, avgqual, minlen, avglen, medianlen, maxlen = self._checkFile(readsfile)
					if "working" in readsfile:
						self._checkPassing(numOfReads, assembled, q30, avgqual)
				else:
					numOfReads, assembled, q30, avgqual, minlen, avglen, medianlen, maxlen = ("FnF", "FnF", "FnF", "FnF", "FnF", "FnF", "FnF", "FnF")
					#outfile.write("\t".join([readsfile, str(numOfReads), str(assembled), str(q30), str(avgqual), str(minlen), str(avglen), str(medianlen), str(maxlen)]) + "\n")
	
				print "\t".join([readsfile, str(numOfReads), str(assembled), str(q30), str(avgqual), str(minlen), str(avglen), str(medianlen), str(maxlen)])
				outfile.write("\t".join([readsfile, str(numOfReads), str(assembled), str(q30), str(avgqual), str(minlen), str(avglen), str(medianlen), str(maxlen)]) + "\n")
		self._setEnvironment()
		return
		
	def _checkFile(self, filename):
		""" """
		qualities = []
		lengths = []
		FQ = FastqStream(filename)
		for fastq in FQ:
			quality = fastq.split('\n')[3].strip()
			qualities.append(average(sanger2phred(quality)))
			lengths.append(len(quality))
		if lengths:
			numOfReads = len(lengths)
			assembled  = round((float(sum(lengths))/ 1000000), 1)
			q30        = ht30(qualities)
			avgqual    = round(average(qualities), 1) 
			minlen     = min(lengths)
			avglen     = int(average(lengths))
			medianlen  = int(median(lengths))
			maxlen     = max(lengths)
			return numOfReads, assembled, q30, avgqual, minlen, avglen, medianlen, maxlen
		else:
			return 0, 0, 0, 0, 0, 0, 0, 0
		
	def _checkPassing(self, numOfReads, assembled, q30, avgqual):
		""" parametrizzare all'esterno """
		criteria1 = (numOfReads >= 40000)
		criteria2 = (assembled >= 50)
		criteria3 = (avgqual >= 22)
		criteria4 = (q30 > 50)
		if criteria1: #and criteria2 and criteria3: # and criteria4:
			self.passing = True
		return
		
	def _setEnvironment(self):
		""" """
		if self.passing == False:
			print "----------------------------------------------------------------------------------------------------------------------------------------------------"
			print "# %s. NOT Passing: %s" % (self.basename, self.passing)
			print "----------------------------------------------------------------------------------------------------------------------------------------------------"
			os.system("mv %s* todonotuse" % self.basename)
			os.system("cp todonotuse/%s_samplereadscheck.tab ." % self.basename) # cambiare qui o spostare dopo
		else:
			print "----------------------------------------------------------------------------------------------------------------------------------------------------"
			print "# %s. Passing: %s" % (self.basename, self.passing)
			print "----------------------------------------------------------------------------------------------------------------------------------------------------"
		return

	def _checkExistence(self):
		"""  """
		for infile in self.readsfiles:
			found = self._which(infile)
			if found:
				self.checked.append(infile)
		return 

	def _which(self, infile):
		if os.path.isfile("./"+infile):
			return True
		return False
		
	def _kurylenko(self):
		""" """
		return
				

### END OF CLASS




if __name__ == "__main__":
	M = ManageSampleReads(sys.argv[1])
	M.manage()
