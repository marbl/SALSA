# SALSA: A tool to scaffold long read assemblies with Hi-C 


This code is used to scaffold your assemblies using Hi-C data. This version implements some improvements in the original SALSA algorithm. If you want to use the old version, it can be found in the `old_salsa` branch. 

To use the latest version, first run the following commands:
```
  cd SALSA
  make
```

To run the code, you will need Python 2.7, [BOOST](http://www.boost.org/) libraries and [Networkx](https://networkx.github.io/).


If you consider using this tool, please cite our publication which describes the methods used for scaffolding.

Ghurye, J., Pop, M., Koren, S., Bickhart, D., & Chin, C. S. (2017). Scaffolding of long read assemblies using long range contact information. BMC genomics, 18(1), 527. [Link](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3879-z)

## How to run the code?

The new version of SALSA has been designed to consider several use cases depending on the input. Some assemblers output assembly graph as well along with the contig sequences. We provide options to use different information provided by the assembly to use for the scaffolding. Here is the what input options look like


### Mapping Reads

To start the scaffolding, first step is to map reads to the assembly. We recommend using BWA or BOWTIE2 aligner to map reads. The read mapping generates a `bam` file. SALSA requires `bed` file as the input. This can be done using the `bamToBed` command from the [Bedtools](http://bedtools.readthedocs.io/en/latest/) package. Also, SALSA requires bed file to be sored by the read name, rather than the alignment coordinates. Once you have bam file, you can run following commands to get the bam
file needed as an input to SALSA.

```
bamToBed -i alignment.bam > alignment.bed
sort -k 4 alignment.bed > tmp && mv tmp alignment.bed
```

### Generating contig lengths file

SALSA requires contig lengths as an input. You can generate this file using this command on your contig sequence file.
```
samtools faidx contigs.fasta
```

This will generate `contigs.fasta.fai` as an output. You can use this file for contig lengths while running SALSA.

### I have contig sequence and the alignment bam file
This is the minimum input you will require Suppose you only have contig sequences generated 
