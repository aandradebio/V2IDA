![alt tag](https://user-images.githubusercontent.com/57667417/84268663-661fda00-aafe-11ea-98c4-abba931a5194.jpg)

# Viral Vaccine genetic Diversity Analyzer 

This script provides an automatized and user-friendly scientific pipeline to perform variant calling and/or quasispecies reconstruction specifically for viral vaccine samples. It was previously used to establish the relationship among genetic diversity, vaccine stability, and the possible reversion to virulence caused by the presence of SNPs and viral quasispecies in vaccine lots from 17DD vaccine against Yellow Fever.

## Pipeline Overview

![alt tag](https://user-images.githubusercontent.com/57667417/84268671-68823400-aafe-11ea-8b58-6fa673230e11.jpg?s=500)

## Setup

### List of Tools Used in this Pipeline

All requirements should be downloaded and installed by the user. 

As default the tools should be in path. As an alternative, the pre compiled files should be in the 00_Exe folder.

[JDK 7](http://jdk7.java.net/)

[BWA](https://github.com/lh3/bwa) v. 

[Samtools](https://github.com/samtools/samtools) v. 

[Picard](https://github.com/broadinstitute/picard) v. 

[GATK](https://github.com/broadinstitute/gatk) v.4

[QuasiRecomb](https://github.com/cbg-ethz/QuasiRecomb) v. 1.2

### Usage

./v2ida.sh ID_NAMES REF I F P
<ul>
<b>ID_NAMES <b> is the name of a .tab file containing the samples names

<b>REF <b> is the reference genome used for alignment (default: fasta format)

<b>I <b> is the inicial nucleotide for quasispecies reconstruction

<b>F <b> is the final nucleotide for quasispecies reconstruction

<b>P <b> is how many parts you would like to divide the genome for quasispecies reconstruction (eg. 1 if you dont want to divide)
<ul>
Example: 

./v2ida.sh sample MN737509 1 10862 5

In this example, the pipeline reads the sample names from the sample.tab file, uses the MN737509.fasta file as reference genome, divides the genome from nucleotide 1 to nucleotide 10.862 in 5 parts. 

Highly customizado. Os filtros são específicos visando amostras de vacinas virais e best practices mas podem ser customizados e associados a outros programas. 

## Credits

This pipeline was developed by [Andrade, AAS](https://github.com/aandradebio) (aandrade@lncc.br) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO), with contributions from [Soares, AER](https://github.com/aersoares81).

ANDRADE, AAS; SOARES, AER; ALMEIDA, LGP; PESTANA, CP; AQUINO, CL; MEDEIROS, MA; VASCONCELOS, ATR. Testing the genomic stability of the Brazilian Yellow Fever vaccine strain using next-generation sequencing data. (IN PRESS)


 


