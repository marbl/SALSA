# SALSA: A tool to scaffold long read assemblies with Hi-C 


This code is used to scaffold your assemblies using Hi-C data. This version implements some improvements in the original SALSA algorithm. If you want to use the old version, it can be found in the `old_salsa` branch. 

To use the latest version, first run the following commands:
```
  cd SALSA
  make
```

To run the code, you will need [BOOST](http://www.boost.org/) libraries and [Networkx](https://networkx.github.io/).

is easy to run on its own as well. You can tweak the parameters in the files to suit your data and run the whole pipeline.

If you consider using this tool, please cite our publication which describes the methods used for scaffolding.

Ghurye, J., Pop, M., Koren, S., Bickhart, D., & Chin, C. S. (2017). Scaffolding of long read assemblies using long range contact information. BMC genomics, 18(1), 527. [Link](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3879-z)
