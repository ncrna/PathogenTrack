# Pathogen-Track
Pathogen-Track is a python-based computational software based on **---** and **---** developped to detect and identify pathogenic microorganisms from single-cell RNA-sequencing (scRNA-seq) raw data. This tool was tested on various scRNA-seq datasets derived from human tumor samples as described in our paper *'Detecting and studying pathogenic microorganisms invasion at the single-cell resolution using Pathogen-Track'*.

Installation
-------------

Before running Pathogen-Track, several dependencies must be installed :

1 . The first step is to install [**UMI-tools**](https://github.com/CGATOxford/UMI-tools). Umi_tools is dependent on python>=3.5, numpy, pandas, scipy, cython, pysam, future, regex and matplotlib, to install it you should start an ssh session and type :

```sh
conda install -c bioconda -c conda-forge umi_tools
```
or
```sh
pip install umi_tools
```


