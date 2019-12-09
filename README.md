#### V2IDA
![alt tag](https://user-images.githubusercontent.com/57667417/69444736-bfc1e200-0d2f-11ea-94f4-592adbb3b0ed.jpeg)

# Viral Vaccine genetic Diversity Analyzer

This script demonstrates using shell scripting to provide an automatized and user-friendly scientific pipeline to perform variant calling and quasispecies reconstruction specifically for viral vaccine samples. 

### List of Tools Used in this Pipeline

[JDK 7](http://jdk7.java.net/)

[BWA](https://github.com/lh3/bwa)

[Samtools](https://github.com/samtools/samtools)

[Picard](https://github.com/broadinstitute/picard)

[GATK](https://github.com/broadinstitute/gatk)

[QuasiRecomb](https://github.com/cbg-ethz/QuasiRecomb)


All requirements should be downloaded and installed by the user. 
As default the tools should be in path.
Or in the 00_Exe folder

### Overview

![alt tag](https://user-images.githubusercontent.com/57667417/69445815-dd904680-0d31-11ea-8885-a2c03c968c92.png)


### Usage

./v2ida.sh samples ref 1 10000 5

where:
samples is the name of a .tab file containing the samples names

ref is the reference genome used for alignment (default: fasta format)

1 is the inicial nucleotide for quasispecies reconstruction

10000 is the final nucleotide for quasispecies reconstruction

5 is how many parts you would like to divide the genome for quasispecies reconstruction


### Credits

This pipeline was developed by Andrade, AAS (aandrade@lncc.br) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO), with contributions from SOARES, AER.

ANDRADE, AAS; SOARES, AER; ALMEIDA, LGP; PESTANA, CP; AQUINO, CL; MEDEIROS, MA; VASCONCELOS, ATR. Testing the genomic stability of the Brazilian Yellow Fever vaccine strain using next-generation sequencing data. (IN PRESS)


 


