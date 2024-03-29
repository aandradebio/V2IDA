<img src="https://user-images.githubusercontent.com/57667417/84274517-599f7f80-ab06-11ea-9ee3-b82e6aa88d75.jpg" width="380">

# Viral Vaccine genetIc Diversity Analyzer 

This script provides an automatized and user-friendly scientific pipeline to perform variant calling and/or quasispecies reconstruction specifically for viral vaccine samples. It was previously used to establish the relationship among genetic diversity, vaccine stability, and the possible reversion to virulence caused by the presence of SNPs and viral quasispecies in vaccine lots from the 17DD strain of Yellow Fever vaccine.

### List of Tools Used in this Pipeline

All requirements should be downloaded and installed by the user. 

As default the tools should be in path. As an alternative, the pre-compiled files should be in the same folder as the V2IDA script.

[JDK 7](http://jdk7.java.net/)

[BWA-MEM](https://github.com/lh3/bwa) v. 0.7

[Samtools](https://github.com/samtools/samtools) v. 1.6

[Picard](https://github.com/broadinstitute/picard) v. 2.21.9

[GATK](https://github.com/broadinstitute/gatk) v.4

[QuasiRecomb](https://github.com/cbg-ethz/QuasiRecomb) v. 1.2

## Pipeline Overview

<img src="https://user-images.githubusercontent.com/57667417/84274511-573d2580-ab06-11ea-9959-ed25f8a5fea2.jpg" width="480">

### Usage
```
./v2ida.sh id ref pair-end initial final parts
```
**id** is the name of a [.tab](https://github.com/aandradebio/V2IDA/blob/master/samples.tab) file containing sample name, r1.fastq file and r2.fastq file in the same line. Other samples should be placed in the next lines. One sample per line; 

**ref** is the reference genome used for alignment (default: fasta format)

**pair-end** or single-end mode

**initial** is the inicial nucleotide for quasispecies reconstruction

**final** is the final nucleotide for quasispecies reconstruction

**parts** is how many parts you would like to divide the genome for quasispecies reconstruction (eg. 1 if you dont want to divide)

**Example:** 
```
./v2ida.sh samples MN737509 pair-end 1 10862 5
```
In this example, the V2IDA pipeline reads the sample names from the samples.tab file, uses the MN737509.fasta file as reference genome, pair-end raw data and divides the genome from nucleotide 1 to nucleotide 10.862 in 5 parts. 

To costumize the SNP hard-filtering criteria, we suggest the reading of [GATK'S Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows). 

Once the V2IDA pipeline analysis is finished, it generates multiple files that comprise metrics and can be opened in any web browser or text editor.

Additionally, we suggest the use of complementary algorithms to perform SNP effect prediction (eg. [SNPeff](https://github.com/pcingola/SnpEff)) and Phylogenetic analysis of reconstructed quasispecies (eg. [Seaview](http://doua.prabi.fr/software/seaview)) from the output files generated by V2IDA pipeline. 


## Credits

This pipeline was developed by [Andrade, AAS](https://github.com/aandradebio) (aandradebio@gmail.com) at the National Laboratory for Scientific Computing - Bioinformatic Laboratory (LABINFO), with contributions from [Soares, AER](https://github.com/aersoares81), Almeida LGP and Vasconcelos, ATR.

Andrade AAS, Soares AER, Paula de Almeida LG, Ciapina LP, Pestana CP, Aquino CL, Medeiros MA, Ribeiro de Vasconcelos AT. Testing the genomic stability of the Brazilian yellow fever vaccine strain using next-generation sequencing data. Interface Focus. 2021 Jun 11;11(4):20200063. doi: 10.1098/rsfs.2020.0063. PMID: 34123353; PMCID: PMC8193464.


 


