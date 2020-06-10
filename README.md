#### V2IDA
![alt tag](https://user-images.githubusercontent.com/57667417/84268663-661fda00-aafe-11ea-98c4-abba931a5194.jpg)

# Viral Vaccine genetic Diversity Analyzer

This script provides an automatized and user-friendly scientific pipeline to perform variant calling and/or quasispecies reconstruction specifically for viral vaccine samples. It was previously was used to establish the relationship among genetic diversity, vaccine stability, and the possible reversion to virulence caused by the presence of SNPs and quasispecies in 17DD vaccine lots.

### List of Tools Used in this Pipeline

[JDK 7](http://jdk7.java.net/)

[BWA](https://github.com/lh3/bwa) v. 

[Samtools](https://github.com/samtools/samtools) v. 

[Picard](https://github.com/broadinstitute/picard) v. 

[GATK](https://github.com/broadinstitute/gatk) v.4

[QuasiRecomb](https://github.com/cbg-ethz/QuasiRecomb) v. 1.2

All requirements should be downloaded and installed by the user. 

As default the tools should be in path.

Or the pre compiled files should be in the 00_Exe folder

### Overview

![alt tag](https://user-images.githubusercontent.com/57667417/84268671-68823400-aafe-11ea-8b58-6fa673230e11.jpg)

### Usage

'''
 ./v2ida.sh ID_NAMES REF I F P
'''

where:
ID_NAMES is the name of a .tab file containing the samples names

REF is the reference genome used for alignment (default: fasta format)

I is the inicial nucleotide for quasispecies reconstruction

F is the final nucleotide for quasispecies reconstruction

P is how many parts you would like to divide the genome for quasispecies reconstruction

Highly customizado. Os filtros são específicos visando amostras de vacinas virais e best practices mas podem ser customizados e associados a outros programas. 

### Credits

This pipeline was developed by [Andrade, AAS](https://github.com/aandradebio) (aandrade@lncc.br) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO), with contributions from [Soares, AER](https://github.com/aersoares81).

ANDRADE, AAS; SOARES, AER; ALMEIDA, LGP; PESTANA, CP; AQUINO, CL; MEDEIROS, MA; VASCONCELOS, ATR. Testing the genomic stability of the Brazilian Yellow Fever vaccine strain using next-generation sequencing data. (IN PRESS)


 


