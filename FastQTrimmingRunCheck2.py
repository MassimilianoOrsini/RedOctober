# -*- coding: utf-8 -*-
"""
This tool trims paired-end FASTQ files on the basis of quality score or left/right position, retaining mate integrity.
Reads without mate after filtering are saved in a separate output file.
"""

import math
#import statistics
import optparse
import sys
import argparse
import time


class PairEndTrimmer:
    
    def __init__(self, input1, input2, maxlengthtrim, lefttrim, righttrim, minqualtrim, avgqualtrim, minlen, trimmed1, trimmed2, trimmedunpaired):# logfile):
        self.input1 = input1
        self.input2 = input2
        self.maxlengthtrim = int(maxlengthtrim)
        self.lefttrim = int(lefttrim)
        self.righttrim = int(righttrim)
        self.minqualtrim = int(minqualtrim)
        self.avgqualtrim = int(avgqualtrim)
        self.minlen = int(minlen)
        self.trimmed1 = trimmed1
        self.trimmed2 = trimmed2
        self.trimmedunpaired = trimmedunpaired
        self.interalaced = input1.split('_R1_001.fastq')[0] +'_working.fastq'
        self.reslinefile = input1.split('_R1_001.fastq')[0]+'_resline.csv'
  
    
    def trimm(self):
        """ """
        total_reads = 0
        discarded_reads = 0
        passing_paired_reads = 0
        passing_unpaired_reads = 0
        forward = open(self.input1)
        reverse = open(self.input2)
        trimmed_forward = open(self.trimmed1, 'w')
        trimmed_reverse = open(self.trimmed2, 'w')
        trimmed_unpaired = open(self.trimmedunpaired, 'w')
        trimmed_interlaced = open(self.interalaced, 'w')
        rawLen = 0
        trimmLen = 0
        rawQual = 0
        trimmQual = 0
        overHangRaw = 0
        overHangTrim = 0
        rawAvgQuals = []
        trimmAvgquals = []
        
        try:
            while True:
                headL = forward.next().rstrip()
                sequL = forward.next().rstrip()
                commL = forward.next().rstrip()
                sangL = forward.next().rstrip()
                qualL = self._sanger2phred(sangL)
                trimmed_sequL, trimmed_qualL = self._trimming(sequL, qualL, self.maxlengthtrim, self.lefttrim, self.righttrim, self.minqualtrim, self.avgqualtrim)
                rawAvgQuals.append(self._average(qualL))
                trimmAvgquals.append(self._average(trimmed_qualL))
                try:
                    headR = reverse.next().rstrip()
                    sequR = reverse.next().rstrip()
                    commR = reverse.next().rstrip()
                    sangR = reverse.next().rstrip()
                except StopIteration:
                    sys.exit('Reverse FASTQ file contain less reads than forward FASTQ file.')
                qualR = self._sanger2phred(sangR)
                trimmed_sequR, trimmed_qualR = self._trimming(sequR, qualR, self.maxlengthtrim, self.lefttrim, self.righttrim, self.minqualtrim, self.avgqualtrim)
                rawAvgQuals.append(self._average(qualR))
                trimmAvgquals.append(self._average(trimmed_qualR))
  
                # Filter by residual length
                if len(trimmed_sequL) >= self.minlen and len(trimmed_sequR) >= self.minlen:
                    trimmed_forward.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + self._phred2sanger(trimmed_qualL) + '\n')
                    trimmed_reverse.write(headR + '\n' + trimmed_sequR + '\n' + commR + '\n' + self._phred2sanger(trimmed_qualR) + '\n')
                    trimmed_interlaced.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + self._phred2sanger(trimmed_qualL) + '\n')
                    trimmed_interlaced.write(headR + '\n' + trimmed_sequR + '\n' + commR + '\n' + self._phred2sanger(trimmed_qualR) + '\n')
                    passing_paired_reads += 2
                elif len(trimmed_sequL) >= self.minlen and len(trimmed_sequR) < self.minlen:
                    trimmed_unpaired.write(headL + '\n' + trimmed_sequL + '\n' + commL + '\n' + self._phred2sanger(trimmed_qualL) + '\n')
                    passing_unpaired_reads += 1
                elif len(trimmed_sequL) < self.minlen and len(trimmed_sequR) >= self.minlen:
                    trimmed_unpaired.write(headR + '\n' + trimmed_sequR + '\n' + commR + '\n' + self._phred2sanger(trimmed_qualR) + '\n')
                    passing_unpaired_reads += 1
                else:
                    discarded_reads += 1
                    
                rawLen = rawLen +len(qualL) + len(qualR)
                #rawQual = rawQual + self._ht30(qualL) + self._ht30(qualR) 

                trimmLen = trimmLen +len(trimmed_sequL) + len(trimmed_sequR)
                #trimmQual = trimmQual + self._ht30(trimmed_qualL) + self._ht30(trimmed_qualR)  
                
                overHangRaw += self._overHangRaw(sequL, sequR)
                overHangTrim += self._overHangTrim(trimmed_sequL, trimmed_sequR)
  
                total_reads += 2
        except StopIteration:
            try:
                reverse.next()
            except StopIteration:
                pass
            else:
                sys.exit('Forward FASTQ file contain less reads than reverse FASTQ file.')

        finally:
            forward.close()
            reverse.close()
            trimmed_forward.close()
            trimmed_reverse.close()
            trimmed_unpaired.close()
            #print "Debug", rawAvgQuals[:50], rawLen
            if float(rawLen) > 0:
				print "DEB", self._ht30(rawAvgQuals)
				q30 = str(int(float(self._ht30(rawAvgQuals))/float(total_reads) *100))
            else:
				q30 = 0
                    
            if float(trimmLen) > 0:
				#print "DEBB", 
				q30Post = str(int(float(self._ht30(trimmAvgquals))/float(passing_paired_reads + passing_unpaired_reads) *100))  
            else:
				q30Post = 0
            
            assembled = str(round((float(rawLen)/ 1000000), 1))
            assembledPost = str(round((float(trimmLen)/ 1000000), 1))        
        
            diffReads = str(int(total_reads) - int(passing_paired_reads))
            diffAssembled = str(float(assembled) - float(assembledPost))
            
            avgqualRaw = str(round((self._median(rawAvgQuals)), 1))  
            avgqualPost = str(round((self._median(trimmAvgquals)), 1))             
            
            resline = [self.input1, self.input2, str(total_reads), assembled, avgqualRaw, q30, str(passing_paired_reads), assembledPost, avgqualPost, q30Post, diffReads, diffAssembled, str(passing_unpaired_reads), str(overHangRaw), str(overHangTrim)]
            with open(self.reslinefile, "wb") as reslinefile:
				reslinefile.write("\t".join(resline) +"\n")

        return

    def _median(self, values):
        """ Median of a list of values """
        theValues = sorted(values)
        if len(theValues) % 2 == 1:
            return int(theValues[(len(theValues)+1)/2-1])
        else:
            try:
                lower = theValues[len(theValues)/2-1]
                upper = theValues[len(theValues)/2]
                return int((float(lower + upper)) / 2)
            except:
                return float('nan')
        
    def _average(self, values):
        """ Arithmetic mean of a list of values """
        return math.fsum(values) / len(values) if len(values) else 0#float('nan')
        
    def _phred2sanger(self, phred_scores):
        """ Convert an array of Phred quality scores (integers in [0, 93]) to a Sanger-encoded quality string"""
        return ''.join([chr(score + 33) for score in phred_scores])
        
    def _sanger2phred(self, sanger_string):
        """ Convert a Sanger-encoded quality string to an array of Phred quality scores"""
        return [ord(ch) - 33 for ch in sanger_string]
        
    def _trimming(self, sequ, qual, maxlengthtrim, lefttrim, righttrim, minqualtrim, avgqualtrim):
        """ Trimming of sequence and quality of a read"""
        # Maximum length trimming
        if maxlengthtrim != -1:
            if len(sequ) > maxlengthtrim:
                sequ = sequ[: maxlengthtrim]
                qual = qual[: maxlengthtrim]
        # Left- and right-side trimming
        if righttrim == 0:
            sequ = sequ[lefttrim :]
            qual = qual[lefttrim :]
        else:
            sequ = sequ[lefttrim : -righttrim]
            qual = qual[lefttrim : -righttrim]
        # Minimum quality right-side trimming
        while len(sequ) and qual[-1] < minqualtrim:
            qual = qual[:-1]
            sequ = sequ[:-1]
        # Average quality right-side trimming
        while len(sequ) and self._average(qual) < avgqualtrim:
            qual = qual[:-1]
            sequ = sequ[:-1]
        return sequ, qual
        
    def _ht30(self, values):
        """ Fraction of values with value higher than 30 """
        tmp = 0
        for element in values:
			try:
				if int(element) >= 30:
					tmp += 1
			except:
				print element
				#sys.exit()
        return tmp
        
    def _overHangRaw(self, read1, read2):
		"""  It Calculates sequence overHang between mates """
		read2 = self._reversecomplement(read2)
		read1 = read1[:-2]
		read2 = read2[2:]
		tail = read1[-5:]
		if read2.startswith(tail):
			return 1
		else:
			return 0
			
    def _overHangTrim(self, read1, read2):
		"""  It Calculates sequence overHang between mates """
		read2 = self._reversecomplement(read2)
		tail = read1[-5:]
		if read2.startswith(tail):
			return 1
		else:
			return 0
					
    def _reversecomplement(self, read):
        """ Return the reverse complement of the dna string."""
        basecomplement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
        nucleotides = list(read)
        nucleotides.reverse()
        nucleotides = [basecomplement[base] for base in nucleotides]
        return ''.join(nucleotides)
        
    def _kurylenko(self):
		""" """
		return
    
    
### END OF CLASS

def main(argv):
    """ main function """
    parser = argparse.ArgumentParser()

    parser.add_argument('-1', dest='input1', action='store', required=True,  help='forward reads file in Sanger FASTQ format')
    parser.add_argument('-2', dest='input2', action='store', required=True, help='reverse reads file in Sanger FASTQ format')
    parser.add_argument('-maxlt', dest='maxlengthtrim', action='store', required=False, default=350, help='maximum length trimming')
    parser.add_argument('-lt', dest='lefttrim', action='store', required=False, default=0, help='left-side trimming')
    parser.add_argument('-rt', dest='righttrim', action='store', required=False, default=0, help='right-side trimming')
    parser.add_argument('-minqt', dest='minqualtrim', action='store', required=False, default=15, help='minimum quality right-side trimming')
    parser.add_argument('-avgqt', dest='avgqualtrim', action='store', required=False, default=20, help='average quality right-side trimming')
    parser.add_argument('-minlf', dest='minlen', action='store', required=False, default=25, help='minimum length filtering')
    parser.add_argument('-t1', dest='trimmed1', action='store', required=True, help='trimmed forward FASTQ file')
    parser.add_argument('-t2', dest='trimmed2', action='store', required=False, help='trimmed reverse FASTQ file')
    parser.add_argument('-tunp', dest='trimmedunpaired', action='store', required=True, help='trimmed unpaired FASTQ file')
    #parser.add_argument('-log', dest='logfile', action='store', required=True, help='log file')    
    
    args = parser.parse_args()

    PET = PairEndTrimmer(args.input1, args.input2, args.maxlengthtrim, args.lefttrim, args.righttrim, args.minqualtrim, args.avgqualtrim, args.minlen, args.trimmed1, args.trimmed2, args.trimmedunpaired)#, args.logfile)
    PET.trimm()
    return

################
if __name__ == "__main__":
    sys.exit(main(sys.argv))
