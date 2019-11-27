#!/bin/bash

###STEP 01 - SETUP AND USER INPUT
## Set Params via Command Line
ID_NAMES=$1
## The name of a file containing the samples names (without the extension)
R=$2
REF=$R.fasta
## Reference genome (default: fasta format)
GSI=$3
## Inicial nucleotide for quasispecies reconstruction
GSF=$4
## Final nucleotide for quasispecies reconstruction
PART=$5
## How many parts you would like to divide the genome (integer number)
WD="$(pwd)"
## Get the current working directory
EX=$WD/00_Exe
#Default: the dependencies should be in path. If not, the executables should be in the 00_Exe folder

RD=$WD/01_Raw_Data
AL=$WD/02_Alignment
VC=$WD/03_Variant_Calling
QR=$WD/04_Quasispecies_Reconstruction
mkdir $AL $VC $QR

## Indexing and creating the dictionary
$EX/bwa/bwa index $REF
$EX/samtools/samtools faidx $REF
java -jar $EX/picard.jar CreateSequenceDictionary R=$REF O=$WD/$R.dict

ntd=$(((GSF-GSI)/PART))
u=$((GSF-ntd))

while read ID;
	do

###STEP 02 - ALIGNMENT AGAINST A REFERENCE GENOME
p=1
f=$GSF
		echo "Starting ID $ID"
		$EX/bwa/bwa mem -R '@RG\tID:group1\tSM:ID\tPL:illumina\tLB:lib1\tPU:unit1' $REF $RD/$ID-R1.fastq $RD/$ID-R2.fastq > $AL/aligned_$ID.sam 
		java -jar $EX/picard.jar SortSam INPUT=$AL/aligned_$ID.sam OUTPUT=$AL/sorted_$ID.bam SORT_ORDER=coordinate 
		java -jar $EX/picard.jar CollectAlignmentSummaryMetrics R=$REF I=$AL/sorted_$ID.bam O=$AL/alignment_metrics_$ID.txt
		java -jar $EX/picard.jar CollectInsertSizeMetrics INPUT=$AL/sorted_$ID.bam OUTPUT=$AL/insert_metrics_$ID.txt HISTOGRAM_FILE=$AL/insert_size_histogram_$ID.pdf
                java -jar $EX/picard.jar MarkDuplicates INPUT=$AL/sorted_$ID.bam OUTPUT=$AL/dedup_$ID.bam METRICS_FILE=$AL/metrics_duplicate_$ID.txt
                java -jar $EX/picard.jar BuildBamIndex INPUT=$AL/dedup_$ID.bam
                $EX/samtools/samtools depth $AL/sorted_$ID.bam > $AL/depth_$ID.txt
		echo "Step 02 - Alignment concluded for $ID"

###STEP 03 - VARIANT CALLING

		java -jar $EX/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $AL/dedup_$ID.bam -o $AL/realignment_targets_$ID.list
                java -jar $EX/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $AL/dedup_$ID.bam -targetIntervals $AL/realignment_targets_$ID.list -o $AL/realigned_reads_$ID.bam
                java -jar $EX/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF -I $AL/realigned_reads_$ID.bam -o $VC/raw_variants_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $VC/raw_variants_$ID.vcf -selectType SNP -o $VC/raw_snps_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $VC/raw_variants_$ID.vcf -selectType INDEL -o $VC/raw_indels_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V $VC/raw_snps_$ID.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $VC/filtered_snps_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V $VC/raw_indels_$ID.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $VC/filtered_indels_$ID.vcf
		java -jar $EX/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $AL/realigned_reads_$ID.bam -knownSites $VC/filtered_snps_$ID.vcf -knownSites $VC/filtered_indels_$ID.vcf -o $VC/recal_data_$ID.table
                java -jar $EX/GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF -I $AL/realigned_reads_$ID.bam -knownSites $VC/filtered_snps_$ID.vcf -knownSites $VC/filtered_indels_$ID.vcf -BQSR $VC/recal_data_$ID.table -o $VC/post_recal_data_$ID.table
                java -jar $EX/GenomeAnalysisTK.jar -T PrintReads -R $REF -I $AL/realigned_reads_$ID.bam -BQSR $VC/recal_data_$ID.table -o $AL/recal_reads_$ID.bam
                java -jar $EX/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $REF -I $AL/recal_reads_$ID.bam -o $VC/raw_variants_recal_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $VC/raw_variants_recal_$ID.vcf -selectType SNP -o $VC/raw_snps_recal_$ID.vcf
		java -jar $EX/GenomeAnalysisTK.jar -T SelectVariants -R $REF -V $VC/raw_variants_recal_$ID.vcf -selectType INDEL -o $VC/raw_indels_recal_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V $VC/raw_snps_recal_$ID.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $VC/filtered_snps_final_$ID.vcf
                java -jar $EX/GenomeAnalysisTK.jar -T VariantFiltration -R $REF -V $VC/raw_indels_recal_$ID.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $VC/filtered_indels_final_$ID.vcf
                echo "Step 03 - Variant Calling concluded for $ID"


###STEP 04 - QUASISPECIES RECONSTRUCTION

		while [ $p -le $PART ]; 
			do
			##$EX/samtools/samtools index $AL/sorted_$ID.bam
			##java -jar $EX/QuasiRecomb.jar -i $AL/sorted_$ID.bam -o $QR/$ID/$ID_$u-$f -r $u-$f
			echo "$u - $f"
			((f=$u))
			u=$((u-ntd))
			((p=$p+1))
			done
p=1
u=$((GSF-ntd))
f=$GSF

		echo "Step 04 - Quasispecies reconstruction concluded for $ID"
	done < $ID_NAMES.tab













