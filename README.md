# RedOctober
RedOctober is a project that has been created with the aim of developing a collection of python/bash scripts to manage Illumina sequencer output; with particular emphasis in microbiology applications


Motivation:
Illumina sequencers write output files in a single folder (assuming a previous bcl2fastq conversion).
Frequently, sequencing runs contain homogeneous samples, let's say as example all bacterial...
The aim of RedOctober is to perform a general workflow to manage all samples taking advantage of multi core resources.
This includes reads trimming, quality checking, hoverang assessment, then a denovo assembling, reports mailing. In addition 
an automatic specie assignment is performed on a basis of a custom database. This will allow, when possible, to
split samples in subfolders which will be subjected to organism-dedicated processing pipelines that will be expored in 
specific sub-projects.


Requirements:
Python 2.7: 
FastQC:
pear: 
abyss:
SPaDEs:
samtools:
Bowtie2:



Usage:





