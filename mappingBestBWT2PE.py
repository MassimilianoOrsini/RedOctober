import os
import sys
import string

def getBestTemplate(bestfile):
	""" """
	infile = open(bestfile, 'r').readlines()
	template = infile[1].split()[0].strip()# + ".genomes.fasta"
	print "# file %s template : %s " % (bestfile, template)
	return template

def mappingBWT2(bestfile, template, pileuping, sensitive): 
	""" """
	prefix = bestfile.rsplit("_spades_scaffolds.best", 1)[0]# +".draft"  
	target = "Genomes/"+template  #.split(".genomes.fasta")[0]#+".target"
	template = "Genomes/"+template + ".genome.fasta"
	bwt2dir = ""
	
	os.system(bwt2dir +"bowtie2 --very-sensitive -a --np 0.2 --rdg 6,4 --rfg 6,4 -x %s -1 %s -2 %s -S %s 2>%s" % (target, prefix+".t1.fq", prefix+".t2.fq", prefix+".sam", prefix+".mapping.log"))
	os.system("samtools faidx %s 2>>%s" % (template, prefix+".mapping.log"))
	os.system("samtools view -bt %s %s > %s 2>>%s" % (template, prefix+".sam", prefix+".bam", prefix+".mapping.log"))
	os.system("samtools sort %s %s 2>>%s" % (prefix+".bam", prefix+".sorted", prefix+".mapping.log"))
	os.system("samtools index %s 2>>%s" % (prefix+".sorted.bam", prefix+".mapping.log"))
	os.system("samtools mpileup -uf %s %s > %s 2>>%s" % (template, prefix+".sorted.bam", prefix+".sorted.bcf", prefix+".mapping.log"))
	os.system("bcftools view -cg %s > %s 2>>%s" % (prefix+".sorted.bcf", prefix+".sorted.vcf", prefix+".mapping.log"))

	os.system("rm -f %s %s %s" % (prefix+".sam", prefix+".bam", prefix+".sorted.bcf"))
	return
		
		
if __name__ == "__main__":
	#readsfile = sys.argv[1]
	#bestfile = sys.argv[2]
	#template = getBestTemplate(bestfile)
	bestfile = sys.argv[1]
	template = getBestTemplate(bestfile)
	#pileuping = sys.argv[3]
	#sensitive = sys.argv[4]
	#mappingBWT2(readsfile, template, pileuping, sensitive)
	mappingBWT2(bestfile, template, "Y", "Y")
