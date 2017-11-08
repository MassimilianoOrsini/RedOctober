#! /bin/bash

clear

# Last Update 21 03 2017
# Gattaca Variant
# USAGE:
# bash RedOctober.sh <runpath> <runname> 
# es:
# bash RedOctober.sh /mnt/data/illumina/fastq_out/160422_NS500787_0030_AH5KYCAFXX run30 

# DEFINING PATHS

pyscriptsdir="/mnt/biowork/scriptsM" 
#
refdir="/mnt/biowork/scriptsM/references" 
speciesdir="/mnt/biowork/scriptsM/REFS" 

flaARefPage="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_campylobacter_seqdef&page=downloadAlleles&locus=CAMP1255"
flaAfile="flaA.tfa"
flaAdbname="FlaADB.fasta"

mlst9RefPage="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_brucella_seqdef&page=downloadAlleles&locus="
mlst9file=".tfa"
mlst9dbname="AllBruc9.mfa" 
mlst9ProfilesPage="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_brucella_seqdef&page=downloadProfiles&scheme_id=1"
mlst21ProfilesPage="http://pubmlst.org/perl/bigsdb/bigsdb.pl?db=pubmlst_brucella_seqdef&page=downloadProfiles&scheme_id=2"


# DO NOT CHANGE ANYTHING BELOW THIS LINE !!!!
########################################################################################################


# SET variables
runpath=$1
runname=$2

numofproc=30
spadesproc=8
pearproc=16

source /mnt/script/applicativi/prll/prll.sh

start=$(date)

########################################################################################################
# DEFINING FUNCTIONS
########################################################################################################


# reads managing

trimfn() { python $pyscriptsdir/FastQTrimmingRunCheck2.py -1 $1 -2 $(basename $1 _R1_001.fastq)"_R2_001.fastq" -lt 20 -rt 2 -minqt 24 -avgqt 28 -minlf 70 -t1 $(basename $1 _R1_001.fastq)"_R1.treads.fastq" -t2 $(basename $1  _R1_001.fastq)"_R2.treads.fastq" -tunp $(basename $1  _R1_001.fastq)"_unp.treads.fastq" ; echo "-lt 20 -rt 2 -minqt 24 -avgqt 28 -minlf 70" > $(basename $1  _R1_001.fastq)"_trimming.cfg"  ; }

fastqcfn() { fastqc $1 >/dev/null 2>&1 ; }

pearingfn () { pear -f $1 -r $(basename $1 _R1.treads.fastq)"_R2.treads.fastq" -o $(basename $1 _R1.treads.fastq)"_m" > $(basename $1 _R1.treads.fastq)"_pear.log" ; }

pearingRawfn () { pear -f $1 -r $(basename $1 _R1_001.fastq)"_R2_001.fastq" -o $(basename $1 _R1_001.fastq)"_M" > $(basename $1 _R1_001.fastq)"_pear.log" ; }

convertingqtoafn() { python $pyscriptsdir/FastQtoA.py -i $1 ; }

checkingreadsfn() { python $pyscriptsdir/ManageSampleReads2.py $(basename $1 _resline.csv) > $(basename $1 _resline.csv)"_checkreads.tmp" ; }

silabreadsfn() { python $pyscriptsdir/silab_checkReads2.py $1 ; }


# mlst - typing

listmlstfn() { srst2 --input_se $1 --mlst_db Listeria_monocytogenes.fasta --mlst_definitions lmonocytogenes.txt --mlst_delimiter '_' --report_new_consensus --output $(basename $1 _working.fastq)"_srst2_lmMLST" ; }

campmlstfn() { srst2 --input_se $1 --mlst_db Campylobacter_jejuni.fasta --mlst_definitions campylobacter.txt --mlst_delimiter '_' --report_new_consensus --output $(basename $1 _working.fastq)"_srst2_cjMLST" ; }

salmmlstfn() { srst2 --input_se $1 --mlst_db Salmonella_enterica.fasta --mlst_definitions salmonella.txt --mlst_delimiter '_' --report_new_consensus --output $(basename $1 _working.fastq)"_srst2_seMLST" ; }

flafn() { srst2 --input_se $1 --mlst_db $flaAdbname --mlst_definitions campyloflaA.txt --mlst_delimiter '_' --report_new_consensus --output $(basename $1 _working.fastq)"_srst2_cjFlaA" ; }

mlst9fn() { srst2 --input_se $1 --mlst_db $mlst9dbname --mlst_definitions brucmlst9.txt --mlst_delimiter '_' --report_new_consensus --output $(basename $1 _working.fastq)"_brucmlst9" > $(basename $1 _working.fastq)"_brucmlst9.log" ; }

# assembling

abyssfn41() { samplename=$(basename $1 _working.fastq); 
	ABYSS -k 41 -c 4 -e 2 -E 4 -o $samplename.contigs $1 >$1.abyss.log ; 
	python $pyscriptsdir/FilterFastaByLen.py  $samplename.contigs 200 ; 
	mv $samplename.contigs.filt.fasta $samplename"_abyss_contigs.fasta" ; }

abyssfn17() { samplename=$(basename $1 _working.fastq); 
	ABYSS -k 17 -c 4 -e 2 -E 4 -o $samplename.contigs $1 >$1.abyss.log ; 
	python $pyscriptsdir/FilterFastaByLen.py  $samplename.contigs 200 ; 
	mv $samplename.contigs.filt.fasta $samplename"_abyss_contigs.fasta" ; }

abyssfnPE() { samplename=$(basename $1 _working.fastq); 
	ABYSS -k 17 -c 4 -e 2 -E 4 -o $samplename.contigs $1 >$1.abyss.log ; 
	python $pyscriptsdir/FilterFastaByLen.py  $samplename.contigs 200 ; 
	mv $samplename.contigs.filt.fasta $samplename"_abyss_contigs.fasta" ; }


AbyssFn()  { samplename=$(basename $1 .fastq)
	for k in 17 21 31  
		do
			ABYSS -k $k -c 4 -e 2 -E 4 -o $samplename"_"$k.fa $1 >> runabyss.log
		done
	cat $samplename"_"??.fa > $samplename"_all.fasta"
	rm $samplename"_"??.fa
	python $pyscriptsdir/FilterFastaByLen.py  $samplename"_all.fasta" 200
	python $pyscriptsdir/nr100.py $samplename"_all.fasta.filt.fasta"
	mv $samplename"_all.fasta.filt.nr.fas" $(basename $samplename "_working")"_abyss_contigs.fasta"
	rm $samplename"_all.fasta.filt.fasta" $samplename"_all.fasta" ; }



spadesSefn() { spades.py --only-assembler -k 21,33,55,77 -t 4 -o $1"_spades" -s $1 >$1.spadesSE.log 2>&1 ; 
	cp $1"_spades"/"scaffolds.fasta" $(basename $1 _working.fastq)".spades.scaffolds" ; 
	rm $1.spadesSE.log ; }

spadesPefn() { spades.py --only-assembler -k 21,33,55,77 -t 4 -o $(basename $1 _R1.treads.fastq)"_spades" -1 $1 -2 $(basename $1 _R1.treads.fastq)"_R2.treads.fastq" -s $(basename $1 _R1.treads.fastq)"_m.merg.fastq" > $1.spadesPE.log 2>&1 ; 
	cp $(basename $1 _R1.treads.fastq)"_spades"/"scaffolds.fasta" $(basename $1 _R1.treads.fastq)".spades.scaffolds"; 
	rm $1.spadesPE.log 
	echo "spades.py --only-assembler -k 21,33,55,77 -t 4  -1 -2 -s " > $1.spadesPE.cfg ; }


filtfastafn() { python $pyscriptsdir/FilterFastaByLen.py  $1 200 ; }


checkdratsfn() { python $pyscriptsdir/CheckDraftLinear.py -d $1 -g 5.2 -m iupac ; }
checkCampyloContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 1.9 -m iupac > $(basename $1 .fasta)".summary" ; }
checkListeriaContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 2.9 -m iupac > $(basename $1 .fasta)".summary" ; }
checkSalmonellaContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 5.2 -m iupac > $(basename $1 .fasta)".summary" ; }
checkBrucellaContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 3.0 -m iupac > $(basename $1 .fasta)".summary" ; }
checkMycoplasmaContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 1.1 -m iupac > $(basename $1 .fasta)".summary" ; }


checkContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 4.0 -m iupac > $(basename $1 .fasta)".summary" ; }
checkRoughContigsfn() { python $pyscriptsdir/CheckContigsLinear.py -c $1 -g 4.0 -m iupac -l "rough" > $(basename $1 .fasta)".summary" ; }

# annotation

renameScafffn() { python $pyscriptsdir/RenameScaffolds.py $1 ; }

prodigalfn() { prodigal -i $1 -f gff -c -q -a $(basename $1 _scaffolds.fasta)"_prodigal_annot.faa" -o $(basename $1 _scaffolds.fasta)"_prodigal_annot.gff" ; }

barrnapfn() { barrnap --quiet $1 > $(basename $1 _scaffolds.fasta)"_barrnap_annot.gff" ; }

ribosomalseqfn() { python $pyscriptsdir/Get16Ssequence.py -s $1 -b $(basename $1 _scaffolds.fasta)"_barrnap_annot.gff" ; }

renameAnnotation() { python $pyscriptsdir/RenameAnnotationFasta.py $(basename $1 _scaffolds.fasta)"_prodigal_annot.faa" ; }

# guessing specie

johnDoefn() { python $pyscriptsdir/FindTemplate.py -i $1 -t $speciesdir/johndoedb -o $(basename $1 _scaffolds.fasta)".species" -x ATGAC -w  > $(basename $1 _scaffolds.fasta)".species.log" ; }
johnDoefnC() { python $pyscriptsdir/FindTemplate.py -i $1 -t $speciesdir/johndoedb -o $(basename $1 _contigs.fasta)".speciesC" -x ATGAC -w  > $(basename $1 _contigs.fasta)".speciesC.log" ; }


findListeriafn() { python $pyscriptsdir/FindTemplate.py -i $1 -t $refdir/Genomi.db -o $(basename $1 _scaffolds.fasta)".best" -x ATGAC -w > $(basename $1 _scaffolds.fasta)".besttempl.log" ; }


# drafting by mapping

mappingBestfn() { python $pyscriptsdir/mappingBestBWT2.py $1 $(basename $1 "_working.fastq")".best" 2>>"assembling00.log" ; }

convertvcf2fqfn() { perl $pyscriptsdir/vcfutils.pl vcf2fq $1 > $(basename $1 .ds_map.sorted.vcf)".cns.fastq" ; }

convertcnsfq2fafn() { python $pyscriptsdir/FastQtoAmultiline.py -q $1 ; }

convertcnsfn() { python $pyscriptsdir/FastQtoAmultiline.py -q $1 ; }

# varie
zippingfn() { tar -czf $1".tar.gz" $1 ; }

btvfiltfastafn() { python $pyscriptsdir/FilterFastaByLen.py  $1 400 ; }

assignsegfn() { python $pyscriptsdir/CheckBtvDenovoAssembly.py -i $1 ; }






########################################################################################################
# RUNNING
########################################################################################################

echo "###############################################################"
echo "# RED OCTOBER PipeLine Gattaca   "
echo "###############################################################"
echo "# Run : " $runname
echo "# Using " $numofproc  " processors"
echo "# Using " $spadesproc  " processors for spades"
echo "# START pipeline :  " $start
echo "###############################################################"


# SET ENVIRONMENT
########################################################################################################
echo "# Set Environment         "
mkdir $runname
cd $runname
echo $runpath > $runname"_sourcedir"
#mkdir tobioinfonas
mkdir todonotuse
mkdir unusedfiles
# Get List dir. 
ls $runpath/*".fastq" | grep -v "Und" > $runname"_listafiles"
cp $runpath/SampleSheet.csv .

# Check Unzipped
ls $runpath/*".fastq.gz" > $runname"_gzippedReads" 2>/dev/null
if [ -s $runname"_gzippedReads" ]; then
	unzipped=$(ls $runpath/*".fastq.gz" | wc -l)
	echo "# ALERT:" $unzipped " UNZIPPED files Found!"
	echo "ciao ecco i file unzipped "$2 | mail -s " Unzipped file list "$2 -a $runname"_gzippedReads" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com
fi

# Set Symlinks
while read A; do ln -s $A; done < $runname"_listafiles"
# mailing START
echo "ciao ecco il sample list corsa "$2 | mail -s " START samples list "$2 -a $runname"_listafiles" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com


# MANAGING READS 
########################################################################################################
now=$(date)
echo "# Checking Quality           " $now
prll -Q -c $numofproc fastqcfn *"R1_001.fastq"  >> $runname"_logfile"
prll -Q -c $numofproc fastqcfn *"R2_001.fastq"  >> $runname"_logfile"
rm *_fastqc.zip                                                        # capire se possibile oppure recuperare il file txt per il trimming dinamico
for file in *R1_001_fastqc.html; do mv $file $(basename $file _001_fastqc.html)"_fastqc.html"; done
for file in *R2_001_fastqc.html; do mv $file $(basename $file _001_fastqc.html)"_fastqc.html"; done

now=$(date)
echo "# Trimming files             " $now                              # unire il trim e il check
prll -Q -c $numofproc trimfn *"R1_001.fastq"  >> $runname"_logfile"

now=$(date)
echo "# Checking Quality Trimmed   " $now                              # ma serve
prll -Q -c $numofproc fastqcfn *".treads.fastq"  >> $runname"_logfile"
rm *_fastqc.zip  

now=$(date)
echo "# Check Run                  " $now                              # unire il trim e il check
python $pyscriptsdir/CheckRun2.py $runname *"_resline.csv" 


now=$(date)
echo "# Pearing files              " $now
prll -Q -c $pearproc pearingfn *"_R1.treads.fastq"  >> $runname"_logfile"
rm *.discarded.fastq
for file in *.unassembled.forward.fastq; do mv $file $(basename $file .unassembled.forward.fastq)".F.fastq"; done
for file in *.unassembled.reverse.fastq; do mv $file $(basename $file .unassembled.reverse.fastq)".R.fastq"; done
for file in *.assembled.fastq; do mv $file $(basename $file .assembled.fastq)".merg.fastq"; done
rm *_pear.log

prll -Q -c $pearproc pearingRawfn *"_R1_001.fastq"  >> $runname"_logfile"
rm *.discarded.fastq
for file in *.unassembled.forward.fastq; do mv $file $(basename $file .unassembled.forward.fastq)".RawF.fastq"; done
for file in *.unassembled.reverse.fastq; do mv $file $(basename $file .unassembled.reverse.fastq)".RawR.fastq"; done
for file in *.assembled.fastq; do mv $file $(basename $file .assembled.fastq)".RawMerg.fastq"; done
rm *_pear.log


now=$(date)
echo "# Checking Reads             " $now
# SILAB check Reads
prll -Q -c $numofproc silabreadsfn *"_resline.csv"
cat *"_SilabReadsRulesCheck.out" | sort -u > $runname"_SilabReadsRulesCheck.tmp"
tail -n 1 $runname"_SilabReadsRulesCheck.tmp" > $runname"_SilabReadsHead"
cat $runname"_SilabReadsRulesCheck.tmp" >> $runname"_SilabReadsHead"
head -n -1 $runname"_SilabReadsHead" > $runname"_SilabReadsRulesCheck.csv"
rm $runname"_SilabReadsRulesCheck.tmp" $runname"_SilabReadsHead"
paste $runname"_runCheck.csv" $runname"_SilabReadsRulesCheck.csv" > $runname"_runCheck2.csv"
mv $runname"_runCheck2.csv" $runname"_runCheck.csv"
# Reads Report
prll -Q -c $numofproc checkingreadsfn *"_resline.csv"  >> $runname"_logfile"
cat *"_checkreads.tmp" >> $runname"_workingReadsReport.csv"
rm *"_checkreads.tmp" $runname"_runCheck.csv" $runname"_SilabReadsRulesCheck.csv"


# mailing il sample reads check
now=$(date)
echo "# Mailing Run Summaries      " $now
echo "ciao ecco il summary delle reads corsa "$2 | mail -s "working reads "$2 -a $2"_workingReadsReport.csv" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com
echo "ciao ecco il summary della corsa "$2 | mail -s "summary "$2 -a $2"_runCheck.csv" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com
#echo "ciao ecco i low qual della corsa "$2 | mail -s "lowqual "$2 -a $2"_lowQualfiles.csv" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com


# reads Q to A in caso di ksnp3
########################################################################################################
now=$(date)
echo "# Converting Reads           " $now
prll -Q -c $numofproc convertingqtoafn *"_working.fastq" >> $runname"_logfile"
prll -Q -c $numofproc zippingfn *".fa" >> $runname"_logfile"
rm *".fa"


# ACCURATE CONTIGS 
########################################################################################################
now=$(date)
echo "# Preparing Accurate Contigs " $now
prll -Q -c $spadesproc spadesPefn *"_R1.treads.fastq" >> $runname"_logfile" 2>>$runname"_logfile"
prll -Q -c $numofproc filtfastafn *".scaffolds"  >> $runname"_logfile"

for file in *.scaffolds.filt.fasta; do mv $file $(basename $file .spades.scaffolds.filt.fasta)"_spades_scaffolds.fasta"; done

rm *".spades.scaffolds"
rm -rf *"_spades"

prll -Q -c $numofproc checkContigsfn *"_scaffolds.fasta" >> $runname"_logfile"
cat *"_scaffolds.summary" | sort -ur >> $runname"_ScaffoldsRestab.csv"
# mailing il sample reads check
#echo "ciao ecco il summary di SPADES della corsa "$2 | mail -s "SPADES Results "$2 -a $2"_scaffoldsrestab.csv" -r m.orsini@izs.it m.orsini.bioinfo@gmail.com
# check 

# GUESS SPECIES 
########################################################################################################
now=$(date)
echo "# Guessing species           " $now
#SPECIES = ['listeria_m', 'campylo_c', 'salmonella_e', 'btv', 'usutu', 'wnv', 'brucella_m', \
#	'brucella_a', 'campylo_j', 'ecoli', 'mycoplasma_m', 'l_innocua', 'l_ivanovii', 'hsapiens', \
#	'ahsv', 'cav', 'cdv', 'cov', 'ehv', 'pprv', 'prv', 'ptero', 'rfv', 'femv', 'tbe'] 
prll -Q -c $numofproc johnDoefn *"_scaffolds.fasta" >> $runname"_logfile"

for file in *.species
  do 
    python $pyscriptsdir/ParseSpeciesFile.py $file  >> $runname"_Speciesfile.tmp"
  done
cat $runname"_Speciesfile.tmp" | sort -u > $runname"_Speciesfile.txt"
rm $runname"_Speciesfile.tmp"

# improve Specie File
python $pyscriptsdir/ImproveSpeciesFile.py $runname"_Speciesfile.txt" $runname"_runCheck.csv" > $runname"_SpeciesCovFile.csv"
# improve Denovo File
python $pyscriptsdir/ImproveSpeciesFile2.py $runname"_Speciesfile.txt" $runname"_ScaffoldsRestab.csv" > $runname"_ScaffoldsRestab2.csv"
rm $runname"_Speciesfile.txt"
mv $runname"_ScaffoldsRestab2.csv" > $runname"_ScaffoldsRestab.csv"
# cambiare questo
rm *_samplereadscheck.tab

# RUNNING SINGLE PIPELINES 
########################################################################################################
#SPECIES = ['listeria_m', 'campylo_c', 'salmonella_e', 'btv', 'usutu', 'wnv', 'brucella_m', \
#	'brucella_a', 'campylo_j', 'ecoli', 'mycoplasma_m', 'l_innocua', 'l_ivanovii', 'hsapiens', \
#	'ahsv', 'cav', 'cdv', 'cov', 'ehv', 'pprv', 'prv', 'ptero', 'rfv', 'femv', 'tbe'] 
	 
# PENDING = ['bos', 'ovis', 'pedv', 'ehdv' ] 



now=$(date)
if [ -d "btv" ]; then
	(cd "btv"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER BTV PipeLine                           "
	#echo "# Target Folder : btv" 
	echo "# Found : " $numofsamples " BTV samples "	
	echo "# nothing to do for btv yet... ;-) "	
	# LAAP
    cd ..) &
fi



now=$(date)
if [ -d "campylo_j" ]; then
	cd "campylo_j"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Campylobacter PipeLine                           "
	#echo "# Target Folder : campylo_j" 
	echo "# Found : " $numofsamples " campylobacter jejuni samples " 
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" &>>"CampyloPipeline.log"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta" &>>"CampyloPipeline.log"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	 &>>"CampyloPipeline.log"
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta" &>>"CampyloPipeline.log"	 
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" &>>"CampyloPipeline.log"

	### MLST ###
	echo "# MLST typing            " $now
	getmlst.py --species "Campylobacter jejuni" &>>"CampyloPipeline.log"
	### FlaA ###
	curl $flaARefPage -o $flaAfile &>>"CampyloPipeline.log"
	# CLAUDIO PIPE LINE 
	### Draft Improving ###
	### Annotation ###
	### Exiting ###
    cd ..
fi

now=$(date)
if [ -d "campylo_c" ]; then
	# iniziare qui il subshell
	cd "campylo_c"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Campylobacter PipeLine                           "
	#echo "# Target Folder  : campylo_c"  
	echo "# Found : " $numofsamples " campylobacter coli samples "
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" &>>"CampyloPipeline.log"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta" &>>"CampyloPipeline.log"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	 &>>"CampyloPipeline.log"
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta" &>>"CampyloPipeline.log"	 
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" &>>"CampyloPipeline.log"
	
	### MLST ###
	echo "# MLST typing            " $now
	getmlst.py --species "Campylobacter jejuni" &>>"CampyloPipeline.log"
	### FlaA ###
	curl $flaARefPage -o $flaAfile &>>"CampyloPipeline.log"
	# CLAUDIO PIPE LINE 
	### Draft Improving ###
	### Exiting ###
    cd ..	
    # terminare qui il subshell
fi


now=$(date)
if [ -d "brucella_m" ]; then
	cd "brucella_m"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Brucella melit. PipeLine                           "
	#echo "# Target Folder  : brucella_m"  
	echo "# Found : " $numofsamples " brucella samples "
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline

	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	
	
	### Mlst-9 ###
	### Draft Improving ###
	### Exiting ###
    cd ..	
fi


now=$(date)
if [ -d "salmonella" ]; then
	cd "salmonella"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Salmonella PipeLine                           "
	#echo "# Target Folder  : salmonella"  
	echo "# Found : " $numofsamples " Salmonella samples "	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	#prll -Q -c $numofproc checkSalmonellaContigsfn *"_scaffolds.fasta" >> $runname"_logfile"
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline

	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"	
	
	
	### MLST ###
	### Draft Improving ###
	### Exiting ###
    cd ..	
fi

now=$(date)
if [ -d "listeria_m" ]; then
	cd "listeria_m"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Listeria PipeLine                           "
	#echo "# Target Folder : listeria_m"
	echo "# Found : " $numofsamples " listeria samples "
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	prll -Q -c $numofproc checkListeriaContigsfn *"_scaffolds.fasta" >> $runname"_logfile"
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"	
	
	### MLST ###   
	now=$(date)
	echo "# MLST typing            " $now
	getmlst.py --species "Listeria monocytogenes" &>>"ListeriaPipeline.log"
	bowtie2-build Listeria_monocytogenes.fasta Listeria_monocytogenes.fasta &>>"ListeriaPipeline.log"  2>&1
	samtools faidx Listeria_monocytogenes.fasta &>>"ListeriaPipeline.log" 2>&1
	prll -Q -c $numofproc listmlstfn *_working.fastq &>>"ListeriaPipeline.log"
	rm *.tfa *.pileup *.bt2 *.fai *.sam 2>>"ListeriaPipeline.log"
	# including clonal complex
	for file in *_results.txt
	do
		python $pyscriptsdir/GetCC.py $file lmonocytogenes.txt  
		mv $file $(basename $file _lmMLST__mlst__Listeria_monocytogenes__results.txt)"_srst2_mlst.txt"
	done	
    ### DRAFT BY MAPPING  ###
	### Draft Improving ###
	### Exiting ###
    cd ..
fi


now=$(date)
if [ -d "hsapiens" ]; then
	(cd "hsapiens"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER H.Sapiens PipeLine                           "
	#echo "# Target Folder : hsapiens"
	echo "# Found : " $numofsamples " human matching "
	
	### Check Data ###
	now=$(date)
	echo "# Removing Host not implemented yet" $now
	# removing host
	### Exiting ###
	cd ..) &
fi


now=$(date)
if [ -d "mycoplasma_m" ]; then
	cd "mycoplasma_m"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Mycoplasma PipeLine                           "
	#echo "# Target Folder : mycoplasma_m"
	echo "# Found : " $numofsamples " mycoplasma_m "
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	prll -Q -c $numofproc checkMycoplasmaContigsfn *"_scaffolds.fasta" >>"Mycoplasma Pipeline"
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" >>"Mycoplasma Pipeline"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta" >>"Mycoplasma Pipeline"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	>>"Mycoplasma Pipeline"
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta" >>"Mycoplasma Pipeline"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta" >>"Mycoplasma Pipeline"
 
    ### DRAFT BY MAPPING  ###
	### Exiting ###
	cd ..
fi

now=$(date)
if [ -d "l_ivanovii" ]; then
	cd "l_ivanovii"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER Listeria PipeLine                           "
	#echo "# Target Folder : l_ivanovii"
	echo "# Found : " $numofsamples " l_ivanovii samples "
	
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	prll -Q -c $numofproc checkListeriaContigsfn *"_scaffolds.fasta" >> $runname"_logfile"
	#prll -Q -c $numofproc checkListeriaContigsfn *"_contigs.fasta" >> $runname"_logfile"
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"	
	
	### MLST ###   
    ### DRAFT BY MAPPING  ###
	### Draft Improving ###
	### Exiting ###
	cd ..
fi



now=$(date)
if [ -d "l_innocua" ]; then
	cd "l_innocua"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED OCTOBER Listeria PipeLine                           "
	#echo "# Target Folder : l_innocua"
	echo "# Found : " $numofsamples " l_innocua samples "
	### Check Data ###
	now=$(date)
	echo "# Checking Data          " $now
	if ! [ -d "notSatisfyingSamples" ]; then
		mkdir notSatisfyingSamples
	fi
	prll -Q -c $numofproc checkListeriaContigsfn *"_scaffolds.fasta" >> $runname"_logfile"
	# inserire un check sugli scaffold magari incrociando dopp con .species e resline
	
	### Annotation ###
	now=$(date)
	echo "# Scaffolds Annotation   " $now
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"
	prll -Q -c $numofproc prodigalfn *"_scaffolds.fasta"
	prll -Q -c $numofproc barrnapfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc ribosomalseqfn *"_scaffolds.fasta"	
	prll -Q -c $numofproc renameScafffn *"_scaffolds.fasta"	
	
	
	### MLST ###   
    ### DRAFT BY MAPPING  ###
	### Draft Improving ###
	### Exiting ###
	cd ..
fi


now=$(date)
if [ -d "manual" ]; then
	(cd "manual"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *.species | wc -l)
	#echo "# RED NOVEMBER manual PipeLine                           "
	#echo "# Target Folder : manual"
	echo "# Found : " $numofsamples " manual samples "
	
	### Exiting ###
	cd ..) &
fi

now=$(date)
if [ -d "todonotuse" ]; then
	(cd "todonotuse"
	cp ../SampleSheet.csv .
	numofsamples=$(ls *_R1_fastqc.html | wc -l)
	#echo "# RED NOVEMBER todonotuse PipeLine                           "
	#echo "# Target Folder : todonotuse"
	echo "# Found : " $numofsamples " todonotuse samples "
	### Exiting ###
	cd ..) &
fi

# CLAUDIO PIPELINE

cd ..
cwd=$(pwd)
echo "# Claudio PipeLine " $cwd"/"$runname"/" $numofproc
python $pyscriptsdir/pipeline_campylobacter.py $cwd"/"$runname"/" $numofproc


wait

# Exiting

rm -rf unusedfiles

# MOVING TO BIOINFONAS
now=$(date)
echo "# Moving To Bioinfonas " $now
bash $pyscriptsdir/2bioinfonas.sh $runname
cp $runname/$runname"_2bioinfonas.exe" .
bash $runname"_2bioinfonas.exe" 2>>FilesToBioinfonas.errors
cat $runname"_2bioinfonas.exe" 2>>FilesToBioinfonas.sh
mv SampleSheet.csv  $runname"_SampleSheet.csv"

end=$(date)
echo "###############################################################"
echo "# END pipeline :" $end
echo "# ..si direbbe che cantino, signore."
echo "###############################################################"
#cd ..
exit 1
