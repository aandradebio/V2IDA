#!/bin/bash

# Developed by Andrade, AAS (aandradebio@gmail.com) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO)
#                                                                           
#  This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.                                         
#                                                                              
#  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.                                
#           
                                                                                                                                           
###STEP 01 - SETUP AND USER INPUT
## Set Paramters via Command Line

if [ "$1" == "-help" ] ; then
    echo
    echo " Usage: ./ `basename $0` id ref pair-end initial final parts"
    echo
    echo " id	is the name of a .tab file containing a sample name, r1.fastq file and r2.fastq file in the same line. Other samples should be in the next lines."
    echo " ref	is the reference genome used for alignment (default: fasta format). Do not input extension, only the name of the file"
    echo " pair-end or single-end"
    echo " initial  is the inicial nucleotide for quasispecies reconstruction"
    echo " final    is the final nucleotide for quasispecies reconstruction"
    echo " parts    is how many parts you would like to divide the genome for quasispecies reconstruction (eg. 1 if you dont want to divide)" 
    echo
    echo " Please visit github.com/aandradebio/V2IDA for more information"
    echo
    exit 0
fi

## The name of a file containing the samples names (without the extension)
ID_NAMES=$1
## Reference genome (default: fasta format)
R=$2
REF=$R.fasta
## Pair-end or Single end alignment mode
MODE="$3"
## Inicial nucleotide for quasispecies reconstruction
GSI=$4
## Final nucleotide for quasispecies reconstruction
GSF=$5
## How many parts you would like to divide the genome (integer number)
PART=$6
## Get the current working directory
WD="$(pwd)"
#Default: the dependencies should be in PATH. If not, you can create an alias (eg. alias gatk='/path/to/gatk-package/gatk') or put the pre-compile executables in the Working directory

if ! [ -x "$(command -v bwa)" ]; then
  echo 'Error: Bwa is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v samtools)" ]; then
  echo 'Error: Samtools is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v ./picard.jar)" ]; then
  echo 'Error: Picard is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v gatk)" ]; then
  echo 'Error: GATK4 is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v ./QuasiRecomb.jar)" ]; then
  echo 'Error: QuasiRecomb is not installed.' >&2
  exit 1
fi

AL=$WD/01_Alignment
VC=$WD/02_Variant_Calling
QR=$WD/03_Quasispecies_Reconstruction
if [ ! -d $AL ]; then mkdir $AL; fi
if [ ! -d $VC ]; then mkdir $VC; fi
if [ ! -d $QR ]; then mkdir $QR; fi

## Indexing and creating the dictionary
INDEX=$WD/$REF.bwt
if [ ! -f $INDEX ]; then bwa index $REF; fi
FAI=$WD/$REF.fai
if [ ! -f $FAI ]; then samtools faidx $REF; fi
DICT=$WD/$R.dict
if [ ! -f $DICT ]; then java -jar picard.jar CreateSequenceDictionary R=$REF O=$R.dict; fi

ntd=$(((GSF-GSI)/PART))
u=$((GSF-ntd))

while read ID FWDID REVID; 
	do 
	echo "Starting $ID"
p=1
f=$GSF

## Alignment 
	if [[ $MODE != "pair-end" ]]; 
	then 
		if [ ! -f $AL/aligned_$ID.sam ]; then bwa mem -R "@RG\tID:group1\tSM:ID\tPL:illumina\tLB:lib1\tPU:unit1" $REF $FWID > $AL/aligned_$ID.sam; fi 
	else 
		if [ ! -f $AL/aligned_$ID.sam ]; then bwa mem -R "@RG\tID:group1\tSM:ID\tPL:illumina\tLB:lib1\tPU:unit1" $REF $FWID $REVID > $AL/aligned_$ID.sam; fi; fi
	if [ ! -f $AL/sorted_$ID.bam ]; then java -jar picard.jar SortSam INPUT=$AL/aligned_$ID.sam OUTPUT=$AL/sorted_$ID.bam SORT_ORDER=coordinate; fi
	if [ ! -f $AL/dedup_$ID.bam ]; then java -jar picard.jar MarkDuplicates INPUT=$AL/sorted_$ID.bam OUTPUT=$AL/dedup_$ID.bam METRICS_FILE=$AL/metrics_duplicate_$ID.txt; fi
	if [ ! -f $AL/dedup_$ID.bai ]; then samtools index $AL/dedup_$ID.bam; fi
	if [ ! -f $AL/depth_$ID.txt ]; then samtools depth $AL/dedup_$ID.bam > $AL/depth_$ID.txt; fi
	echo "Mapping step concluded for $ID"

## SNP Calling

	if [ ! -f $VC/raw_variants_$ID.vcf ]; then gatk HaplotypeCaller -R $REF -I $AL/dedup_$ID.bam -O $VC/raw_variants_$ID.vcf; fi
	if [ ! -f $VC/raw_SNPs_$ID.vcf ]; then gatk SelectVariants -R $REF -V $VC/raw_variants_$ID.vcf --select-type SNP -O $VC/raw_SNPs_$ID.vcf; fi
	if [ ! -f $VC/filtered_SNPs_$ID.vcf ]; then gatk VariantFiltration -R $REF -V $VC/raw_SNPs_$ID.vcf -O $VC/filtered_SNPs_$ID.vcf --filter-expression "QD < 2.0 && FS > 60.0 && MQ < 40.0 && SOR > 4.0" --filter-name "filter_SNPs"; fi
	if [ ! -f $VC/bqsr_snps_$ID.vcf ]; then gatk SelectVariants --exclude-filtered -V $VC/filtered_SNPs_$ID.vcf -O $VC/bqsr_snps_$ID.vcf; fi
	if [ ! -f $VC/raw_indels_$ID.vcf ]; then gatk SelectVariants -R $REF -V $VC/raw_variants_$ID.vcf --select-type INDEL -O $VC/raw_indels_$ID.vcf; fi
	if [ ! -f $VC/filtered_indels_$ID.vcf ]; then gatk VariantFiltration -R $REF -V $VC/raw_indels_$ID.vcf -O $VC/filtered_indels_$ID.vcf --filter-expression "QD < 0.2 && FS > 200.0 && SOR > 10.0" --filter-name "filter_indels"; fi
	if [ ! -f $VC/bqsr_indels_$ID.vcf ]; then gatk SelectVariants --exclude-filtered -V $VC/filtered_indels_$ID.vcf -O $VC/bqsr_indels_$ID.vcf; fi
	if [ ! -f $VC/recal_data_$ID.table ]; then gatk BaseRecalibrator -R $REF -I $AL/sorted_$ID.bam --known-sites $VC/bqsr_snps_$ID.vcf --known-sites $VC/bqsr_indels_$ID.vcf -O $VC/recal_data_$ID.table; fi
	if [ ! -f $VC/recal_reads_$ID.bam ]; then gatk ApplyBQSR -R $REF -I $AL/sorted_$ID.bam --bqsr $VC/recal_data_$ID.table -O $VC/recal_reads_$ID.bam; fi
	if [ ! -f $VC/raw_variants_recal_$ID.vcf ]; then gatk HaplotypeCaller -R $REF -I $VC/recal_reads_$ID.bam -O $VC/raw_variants_recal_$ID.vcf; fi
	if [ ! -f $VC/recal_SNPs_$ID.vcf ]; then gatk SelectVariants -R $REF -V $VC/raw_variants_recal_$ID.vcf --select-type SNP -O $VC/recal_SNPs_$ID.vcf; fi
	if [ ! -f $VC/final_SNPs_$ID.vcf ]; then gatk VariantFiltration -R $REF -V $VC/recal_SNPs_$ID.vcf -O $VC/final_SNPs_$ID.vcf --filter-expression "QD < 2.0 && FS > 60.0 && MQ < 40.0 && SOR > 4.0" --filter-name "filter_SNPs"; fi
	if [ ! -f $VC/recal_indels_$ID.vcf ]; then gatk SelectVariants -R $REF -V $VC/raw_variants_recal_$ID.vcf --select-type INDEL -O $VC/recal_indels_$ID.vcf; fi
	if [ ! -f $VC/final_indels_$ID.vcf ]; then gatk VariantFiltration -R $REF -V $VC/recal_indels_$ID.vcf -O $VC/final_indels_$ID.vcf --filter-expression "QD < 0.2 && FS > 200.0 && SOR > 10.0" --filter-name "filter_indels"; fi
	echo "SNP calling step concluded for $ID"

## Quasispecies reconstruction

	while [ $p -le $PART ]; 
		do
		if [ ! -f $QR/$ID/$u-$f/quasispecies.fasta ]; then java -jar QuasiRecomb.jar -i $AL/dedup_$ID.bam -o $QR/$ID/$u-$f -noGaps -r $u-$f; fi
		((f=$u))
		u=$((u-ntd))
		((p=$p+1))
		done
	p=1
	u=$((GSF-ntd))
	f=$GSF
		echo "Quasispecies reconstruction step concluded for $ID"
done < $ID_NAMES.tab
