import os
import string
import sys

infilename = sys.argv[1]
tmpfilename = infilename+'.tmp'
logfilename = infilename+'.stat'
samplename = infilename.split('_spades_')[0] # attenzione qui

contigcounter = 1

with open(infilename) as infile, open(tmpfilename, 'wb') as tmpfile, open(logfilename, 'wb') as statfile:
	for fasta in infile.read().split(">")[1:]:
		title = fasta.split("\n", 1)[0]
		sequence = fasta.split("\n", 1)[1]
		tmpfile.write(">"+samplename+'_contig'+str(contigcounter)+'\n'+sequence.strip()+'\n')
		statfile.write(samplename+'_contig'+str(contigcounter)+'\t'+title.strip()+'\n')
		contigcounter += 1
	    
os.system("mv %s %s" % (tmpfilename, infilename))
        
