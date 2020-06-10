![alt tag](https://user-images.githubusercontent.com/57667417/84274517-599f7f80-ab06-11ea-9ee3-b82e6aa88d75.jpg)

# Viral Vaccine genetic Diversity Analyzer 

This script provides an automatized and user-friendly scientific pipeline to perform variant calling and/or quasispecies reconstruction specifically for viral vaccine samples. It was previously used to establish the relationship among genetic diversity, vaccine stability, and the possible reversion to virulence caused by the presence of SNPs and viral quasispecies in vaccine lots from 17DD vaccine against Yellow Fever.

### List of Tools Used in this Pipeline

All requirements should be downloaded and installed by the user. 

As default the tools should be in path. As an alternative, the pre-compiled files should be in the same folder as the script.

[JDK 7](http://jdk7.java.net/)

[BWA](https://github.com/lh3/bwa) v. 0.7.17

[Samtools](https://github.com/samtools/samtools) v. 

[Picard](https://github.com/broadinstitute/picard) v. 

[GATK](https://github.com/broadinstitute/gatk) v.4

[QuasiRecomb](https://github.com/cbg-ethz/QuasiRecomb) v. 1.2

## Pipeline Overview

<img src="https://user-images.githubusercontent.com/57667417/84274511-573d2580-ab06-11ea-9959-ed25f8a5fea2.jpg" width="480">

### Input data

This pipeline requires short-read data and a reference consensus genome as input. 
To facilitate multi-sample usage, all raw data should be named as follow: ```sample_name-R1.fastq and sample_name-R2.fastq``` for pair-end reads or ```sample_name.fastq``` for single-end reads. 

### Usage
```
./v2ida.sh ID_NAMES REF I F P
```
**ID_NAMES** is the name of a [.tab](https://github.com/aandradebio/V2IDA/blob/master/samples.tab) file containing the sequencing mode (pair-end or single-end) followed by samples names

**REF** is the reference genome used for alignment (default: fasta format)

**I** is the inicial nucleotide for quasispecies reconstruction

**F** is the final nucleotide for quasispecies reconstruction

**P** is how many parts you would like to divide the genome for quasispecies reconstruction (eg. 1 if you dont want to divide)

**Example:** 
```
./v2ida.sh samples MN737509 1 10862 5
```
In this example, the V2IDA pipeline reads the sample names from the samples.tab file, uses the MN737509.fasta file as reference genome, divides the genome from nucleotide 1 to nucleotide 10.862 in 5 parts. 

To costumize the SNP hard-filtering criteria, we suggest the reading of [GATK'S Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). 

Once the V2IDA pipeline analysis is finished, it generates multiple files that comprise the general statistics and can be opened in any web browser or text editor.

Additionally, we suggest the usage of complementary algorithms to perform SNP effect prediction (eg. [SNPeff](https://github.com/pcingola/SnpEff)) and Phylogenetic analysis of reconstructed quasispecies (eg. [Seaview](http://doua.prabi.fr/software/seaview)). 


## Credits

This pipeline was developed by [Andrade, AAS](https://github.com/aandradebio) (aandrade@lncc.br) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO), with contributions from [Soares, AER](https://github.com/aersoares81).

ANDRADE, AAS; SOARES, AER; ALMEIDA, LGP; PESTANA, CP; AQUINO, CL; MEDEIROS, MA; VASCONCELOS, ATR. Testing the genomic stability of the Brazilian Yellow Fever vaccine strain using next-generation sequencing data. (IN PRESS)


 


