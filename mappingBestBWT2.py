import os
import sys
import string

def getBestTemplate(bestfile):
	""" """
	infile = open(bestfile, 'r').readlines()
	template = "NC_" + infile[1].split("NC_")[1].split(".")[0].split(" ")[0] + ".fna"
	print "# file %s.best template : %s " % (bestfile, template)
	return template

def mappingBWT2(readsfile, template, pileuping, sensitive):
	""" """
	prefix = readsfile.rsplit("_working.fq", 1)[0] +"_draft"
	target = "/mnt/biowork/scriptsM/references/"+template.split(".fna")[0]+".target"
	template = "/mnt/biowork/scriptsM/references/"+template
	
	os.system(bwt2dir +"bowtie2 --very-sensitive -a --np 0.2 --rdg 6,4 --rfg 6,4 -a -x %s -U %s -S %s 2>%s" % (target, readsfile, prefix+".sam", "mapping.log"))
	os.system("samtools faidx %s 2>>%s" % (template, "mapping.log"))
	os.system("samtools view -bt %s %s > %s 2>>%s" % (template, prefix+".sam", prefix+".bam", "mapping.log"))
	os.system("samtools sort %s %s 2>>%s" % (prefix+".bam", prefix+".sorted", "mapping.log"))
	os.system("samtools index %s 2>>%s" % (prefix+".sorted.bam", "mapping.log"))
	os.system("samtools mpileup -uf %s %s > %s 2>>%s" % (template, prefix+".sorted.bam", prefix+".sorted.bcf", "mapping.log"))
	os.system("bcftools view -cg %s > %s 2>>%s" % (prefix+".sorted.bcf", prefix+".sorted.vcf", "mapping.log"))

	os.system("rm -f %s %s %s" % (prefix+".sam", prefix+".bam", prefix+".sorted.bcf"))
	return
		
		
if __name__ == "__main__":
	 
	readsfile = sys.argv[1]
	bestfile = sys.argv[2]
	print "debug", readsfile, bestfile
	template = getBestTemplate(bestfile)
	pileuping = sys.argv[3]
	sensitive = sys.argv[4]
	mappingBWT2(readsfile, template, pileuping, sensitive)
