# BiAlign

Bialignment of RNAs and proteins. This type of alignment supports
evolutionary 'shift' events that allow some incongruence between sequence and 
structure evolution.

The first version of this tool has been described in 

Waldl M., Will S., Wolfinger M.T., Hofacker I.L., Stadler P.F. (2020)
Bi-alignments as Models of Incongruent Evolution of RNA Sequence and
Secondary Structure. In: Cazzaniga P., Besozzi D., Merelli I., Manzoni L.
(eds) Computational Intelligence Methods for Bioinformatics and
Biostatistics. CIBB 2019. Lecture Notes in Computer Science, vol 12313.
Springer, Cham. https://doi.org/10.1007/978-3-030-63061-4_15



## Setup and usage

The software requires compilation by Cython. For this purpose Cython and
Python need to be installed and the tool has to be compiled by

```
python setup.py build_ext --inplace
```

Moreover, for aligning RNAs, the tool requires the Vienna RNA package with Python bindings. We recommend to use the tool under Linux or MacOS and install the prerequisites via conda / bioconda.

### Example calls

#### Bi-alignment of two RNAs with given structure strings
```
./bialign.py UGUAAACAUCCUCGACUGGAAGCUGUGAAGCCACAAAUGGGCUUUCAGUCGGAUGUUUGCA UGUAAACAUCCUACACUCAGCUGUCAUACAUGCGUUGGCUGGGAUGUGGAUGUUUACG --strA "(((((((((((..((((((((((((((....))).....))))))))))))))))))))))" --strB "(((((((((((.(((((((((((............)))))))).))))))))))))))"
```

#### Bi-Alignment of two proteins.
Input is read from files as written by the secondary structure prediction
web server CFSSP (Kumar et al, 2013; http://www.biogem.org/tool/chou-fasman).

```
./bialign.py --filein Examples/DNAPolymerase1_Ecoli_AS-Struktur.txt Examples/DNAPolymerase1_Xanthomonas_AS-Struktur.txt --max_shift 1 --type Protein
```

For further command line parameters that configure modes and
alignment parameters, please refer to the help output of the tool as obtained by
```
./bialign.py --help
```
