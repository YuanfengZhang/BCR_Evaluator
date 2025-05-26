# BCREval
>[!IMPORTANT]
> This is a customized version of BCREval, with parallelization and acceleration powered by C++.
> The original version can be found [here](https://github.com/hqyone/BCR_Evaluator)
> The following readme info is also changed.
BCREval is used to evaluate base conversion rate in WGBS/RRBS/EM-seq experiments. It will support TAPS soon.

Pull requests are welcome.

## Todo
1. Add support for TAPS
2. Distribute this as a pypi package.

## References
Zhou J, Zhao M, Sun Z, Wu F, Liu Y, Liu X, He Z, He Q, He Q. BCREval: a computational method to estimate the bisulfite conversion ratio in WGBS. BMC Bioinformatics. 2020 Jan 31;21(1):38. doi: 10.1186/s12859-019-3334-z. PMID: 32005131; PMCID: PMC6995172.

## Installation
1. create a conda env
```bash
mamba create -n pybind -c conda-forge \
  python=3.13 scipy pybind11 zstandard
```

or
```bash
uv python install 3.13 && cd BCR_Evaluator && \
   uv python pin 3.13 && uv init && uv add pybind11
```

2. compile the c++ module using g++
```bash
g++ -O3 -shared -std=c++11 -fPIC $(python3-config --includes) \
     telomere_parser.cpp \
     -o telomere_parser$(python3-config --extension-suffix)
```
or use modify the cpp file (PyModuleDef -> PYBIND11_MODULE) and compile it with pybind11.
```bash
c++ -O3 -shared -std=c++11 -fPIC $(python3-config --includes) \
     telomere_parser.cpp \
     -o telomere_parser$(python3-config --extension-suffix)
```

## Usage information
Usage:
```bash
python BCREval.py \
  -i {input} \  # can be .fq, .fq.gz, .fq.zst
  -o {output} \  # the parent dir will be created if it does not exist
  -t {threads} \  # optional, default is 1
  -c {chunk_size} \  # optional, default is 1000000
  -v {verbose}  # optional, debug / info / warning / error / critical
  2 > {error.log}  # The stderr could be very long, prepare.
```

## Output format
> [!NOTE]
> The output format is concise and very differnt from the original version.
The output is a tsv file containing the following columns:
+ **input_fname**: The file name of your input file.
+ **strand**: C / G
+ **total_reads**: The total number of reads in the FASTQ file.
+ **telomeric_reads**: The number of reads containing telomeric repeats.
+ **completely_converted_unit_counts**: The number of completely converted telomeric units.
+ **total_n3_reads**: The total number of reads which contain N3 units that has three methylated cytosines.
+ **strand_seq_counts**: The counter of telomeric units of C-strand derived telomeric reads.
+ **r1**: The unconverted/methylated ratio of the first C.
+ **r2**: The unconverted/methylated ratio of the second C.
+ **r3**: The unconverted/methylated ratio of the third C.
+ **<code style="color : red">optimized_conversion_rate</code>**: base conversion rate calculated and optimized by [L-BFGS-B](https://en.wikipedia.org/wiki/Limited-memory_BFGS) which may be more accurate than the original method.
+ **<code style="color : red">estimated_conversion_rate</code>**: 1 - (r1 + r2 + r3) / 3. This is the approach authors use.
+ **<code style="color : red">conversion_rate_r1r2</code>**: r2 / (r1 + r2). Ignore r3 for robustness. Please read the paper mentioned above for details.
+ **<code style="color : red">conversion_rate_r3</code>**: 1 - r3 ** (1 / 3). Only use r3.
