# benchee

Benchee is a SV benchmarking tool. 

## Installation

First you need to clone the repository:

```
git clone https://github.com/Joshtron/benchee
```

Please make sure to install pybedtools. pip and conda commands are listed below:

```
pip install pybedtools
```

```
conda install -c bioconda pybedtools
```

After successfully installing pybedtools, call the setup.py script which will enable you to use the benchee command line 
tool.

```
python /path/to/benchee_project/setup.py develop
```

## Command line tool

```
benshee --query query.vcf --truth truth.vcf benchmark
```
