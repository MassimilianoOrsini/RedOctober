import os
import string
import sys
import csv

genomesizes = {'listeria_m': 2.9, 'campylo_c': 1.6, 'salmonella_e': 5.2, \
				'btv': 0.05, 'pedv': 0.05, 'usutu': 0.05, 'wnv': 0.05, \
				'brucella_m': 3.0, 'brucella_a': 3.0, 'campylo_j': 1.6, \
				'ecoli': 5.0, 'mycoplasma_m': 0.8, 'l_innocua': 2.9, \
				'l_ivanovii': 2.9, 'ahsv': 0.05, 'cav': 0.05, 'cdv': 0.05, \
				'cov': 0.05, 'ehv': 0.05, 'pprv': 0.05, 'prv': 0.05, \
				'ptero': 0.05, 'rfv': 0.05, 'femv': 0.05, 'tbe': 0.05}



speciesfilename = sys.argv[1] # $runname"_Speciesfile.txt" 
runfilename = sys.argv[2] # $2"_runCheck.csv"

# acquire Called bases
called = {}
with open(runfilename) as checkfile:
	checkreader = csv.reader(checkfile, delimiter='\t')
	for row in checkreader:
		if row[0].startswith('#'):
			continue
		else:
			samplebasename = row[0].split("_R1_001.fastq")[0]
			called[samplebasename] = float(row[7])
	
# add values
with open(speciesfilename) as speciesfile:
	speciesreader = csv.reader(speciesfile, delimiter='\t')
	for row in speciesreader:
		#print row
		if row[0].startswith('# filename'):
			headline = row + ["Coverage"]
			print "\t".join(headline)
		else:
			samplebasename = row[0].split(':')[1].strip()
			recognizedas = row[1].strip()
			if samplebasename in called and recognizedas in genomesizes:
				coverage = str(round((called[samplebasename]/genomesizes[recognizedas]), 1))
				resline = [samplebasename] + row[1:] + [coverage]
				print "\t".join(resline)
			elif samplebasename in called and recognizedas not in genomesizes:
				resline = [samplebasename] + row[1:] + ["not_applicable"]
				print "\t".join(resline)
			else:
				resline = [samplebasename] + row[1:] + ["error_in_samplename"]
				print "\t".join(resline)
			
