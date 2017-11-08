# -*- coding: utf-8 -*-
"""
Check bacterial draft
"""

import math
import optparse
import os
import re

def fasta_size(filename):
    """ Determine the length of a FASTA file containing only 1 sequence"""
    tmp = 0
    with open(filename) as fasta:
        fasta.next() # skip header
        try:
            while True:
                row = fasta.next().strip()
                tmp += len(row)
        except StopIteration:
            return tmp


def average(values):
    """ Arithmetic mean of a list of values """
    return math.fsum(values) / len(values) if len(values) else float('nan')


def median(values):
    """ Median of a list of values """
    if len(values) == 0:
        return float('nan')
    sorted_values = sorted(values)
    if len(sorted_values) % 2 == 1:
        return sorted_values[(len(sorted_values)-1) / 2]
    else:
        lower = sorted_values[len(sorted_values)/2 - 1]
        upper = sorted_values[len(sorted_values)/2]
        return float(lower + upper) / 2


def N50(values):
    """ N50 of a list of contig lengths, i.e. the length for which the set of all contigs of that length or longer contains at least half of the total of the lengths of the contigs, and for which the set of all contigs of that length or shorter contains at least half of the total of the lengths of the contigs. When 2 values meets both these criteria then the N50 is their average. """
    if len(values) == 0:
        return float('nan')
    sorted_values = sorted(values)
    assembled = math.fsum(values)
    index = 0
    tmp50 = sorted_values[index]
    while tmp50 < assembled / 2:
        index += 1
        tmp50 += sorted_values[index]
    return sorted_values[index] if tmp50 > assembled / 2 else float(sorted_values[index] + sorted_values[index + 1]) / 2


def NG50(values, genomesize):
    """ NG50 of a list of contig lengths, i.e. the length for which the set of all contigs of that length or longer contains at least half of the given genome size. """
    if len(values) == 0:
        return float('nan')
    rev_sorted_values = sorted(values, reverse=True)
    index = 0
    tmpG50 = rev_sorted_values[index]
    try:
        while tmpG50 < float(genomesize) / 2:
            index += 1
            tmpG50 += rev_sorted_values[index]
        return rev_sorted_values[index]
    except IndexError:
        return 0


def __main__():
    """ main function """
    parser = optparse.OptionParser()
    parser.add_option('-d', dest='draft', help='draft genome')
    parser.add_option('-g', dest='reference', help='reference genome size')
    parser.add_option('-m', dest='mode', help='contig recognition mode iupac|acgt ')
    (options, args) = parser.parse_args()
    if len(args) > 0:
        parser.error('Wrong number of arguments')
    
    draft = options.draft
    reference = options.reference
    mode = options.mode
    
    genomesize = float(reference) * 1000000 #fasta_size(reference)
    fasta = ""
    with open(draft, 'r') as infile:
        infile.next() # skip header
        try:
            while True:
                row = infile.next().strip().lower()
                fasta += row
        except StopIteration:
            fasta.replace("\n", "")
    
    ns = fasta.count("n")
    nsperc =  float(ns) / len(fasta) * 100
    gaps = re.compile('n+')
    gapsiterator = gaps.findall(fasta)
    if gapsiterator != []:
        gapslength = [len(gap) for gap in gapsiterator]
    else:
        gapslength = [0]
    if mode == "iupac":
        contig = re.compile('[acgtrymksbdhvw]+')
    else:
        contig = re.compile('[acgt]+')
      
    contigiterator = contig.findall(fasta)
    contigslength = [len(contig) for contig in contigiterator]
    headline = ["# Sample", "Len.", "Assembled", "AssembledPerc.", "Contigs"]
    resline = ["# "+draft.split(".")[0], str(len(fasta)), str(sum(contigslength)), str(100 - round((nsperc), 1)), str(len(contigiterator)) ]
    if len(contigiterator):
		headline += ["Avg.", "Median", "Min.", "Max.", "Lt200", "N50", "Gaps", "GapsAvg", "GapsMed", "GapsMax", "GapsSt20"]
		contigs_over200 = sum([el >= 200 for el in contigslength])
		gaps_below5 = sum([el <= 20 for el in gapslength])
		resline.append(str(round((average(contigslength)),1)))
		resline.append(str(median(contigslength)))
		resline.append(str(min(contigslength)))
		resline.append(str(max(contigslength)))
		resline.append(str(contigs_over200))
		resline.append(str(N50(contigslength)))
		resline.append(str(len(gapsiterator)))
		resline.append(str( round((average(gapslength)),1)))
		resline.append(str(median(gapslength)))
		resline.append(str(max(gapslength)))
		resline.append(str(gaps_below5))
		print "\t".join(headline)
		print "\t".join(resline)


if __name__ == "__main__":
    __main__()
