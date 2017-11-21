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

```
python run_pipeline.py -h
usage: run_pipeline.py [-h] -a ASSEMBLY -l LENGTH -b BED [-o OUTPUT]
                       [-c CUTOFF] [-g GFA] [-u UNITIGS] [-t TENX] [-e ENZYME]
                       [-i ITER] [-x DUP] [-s EXP] [-m CLEAN]

SALSA Iterative Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        Path to initial assembly
  -l LENGTH, --length LENGTH
                        Length of contigs at start
  -b BED, --bed BED     Bed file of alignments sorted by read names
  -o OUTPUT, --output OUTPUT
                        Output directory to put results
  -c CUTOFF, --cutoff CUTOFF
                        Minimum contig length to scaffold, default=1000
  -g GFA, --gfa GFA     GFA file for assembly
  -u UNITIGS, --unitigs UNITIGS
                        The tiling of unitigs to contigs in bed format
  -t TENX, --tenx TENX  10x links tab separated file, sorted by last columnls
  -e ENZYME, --enzyme ENZYME
                        Restriction Enzyme used for experiment
  -i ITER, --iter ITER  Number of iterations to run, default = 3
  -x DUP, --dup DUP     File containing duplicated contig information
  -s EXP, --exp EXP     Expected Genome size of the assembled genome
  -m CLEAN, --clean CLEAN
                        Set this option to "yes" if you want to find
                        misassemblies in input assembly
```

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

This will generate `contigs.fasta.fai` as an input for `-l` option. You can use this file for contig lengths while running SALSA.

### Restriction Enzyme input

Hi-C experiments can use different restriction enzymes. We use the enzyme frequency in contigs to normalize the Hi-C interaction frequency. You will need to specify which enzyme was used for Hi-C experiment while running SALSA in `-e` option. If multiple enzymes were used, they can specified by separating with comma without space, like `-e GATC,AAGCTT`.


### 1) I have contig sequences and the alignment bam file
This is the minimum input you will require Suppose you only have contig sequences generated. Once you prepare the bed file as described above, the code can be run as follows:
```
python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds 
```

### 2) I have contig sequences and the alignment bam file but also want to use Hi-C data to correct input assembly errors

We also implemented a method in SALSA that can correct some of the errors in the assembly with Hi-C data. To use this method, you need to run following
```
python run_pipeline.py -a contigs.fasta -l contigs.fasta.fai -b alignment.bed -e {Your Enzyme} -o scaffolds -m yes
```

If you want to know what were the locations in the contigs where SALSA found errors, you can look at the `input_breaks` file in the output directory.

