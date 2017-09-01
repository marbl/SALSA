# SALSA: A tool to scaffold long read assemblies with Hi-C 
## (This work is sponsored by Pacific Biosciences, Menlo Park)

This code is used to scaffold your assemblies using Hi-C data. To use the code first run the following command:
```
  make
```

To run the code, you will need [Samtools](http://samtools.sourceforge.net), [Bedtools](http://bedtools.readthedocs.io/en/latest/), [pbcore](https://github.com/PacificBiosciences/pbcore) and [Networkx](https://networkx.github.io/).

After this, there will be two binary files generated, one for ```break_contigs.cpp``` and one for ```triangle_plot.cpp```. Now you can run the code as follows:

The primary file to run the pipeline is run.py which has following options

```
python run.py -h
usage: run.py [-h] -a ASSEMBLY -m MAPPING -d DIR [-b MISSASSEMBLY]

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        assembled contigs
  -m MAPPING, --mapping MAPPING
                        mapping of read to contigs in bam format
  -d DIR, --dir DIR     output directory for results
  -b MISSASSEMBLY, --missassembly MISSASSEMBLY
                        add flag to find and break misassemblies from the
                        contigs
```

The fasta file containing final scaffolds will be generated in your output directory as ```scaffold.fasta```. Each step in the pipeline is easy to run on its own as well. You can tweak the parameters in the files to suit your data and run the whole pipeline.

If you consider using this tool, please cite our publication which describes the methods used for scaffolding.

Ghurye, J., Pop, M., Koren, S., Bickhart, D., & Chin, C. S. (2017). Scaffolding of long read assemblies using long range contact information. BMC genomics, 18(1), 527. [Link](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3879-z)


## Results

We ran this tool on the assembly of NA12878 with N50 = 1.55 Mb and used 725 million hi-c reads for scaffolding. Our final scaffold had N50 = 80 Mb. Here are the dotplots for the chromosomes

![Alt text](chr.png?raw=true "Optional Title")

