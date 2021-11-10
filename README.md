# BiAlign - Bialignment of RNAs and proteins

This type of alignment supports evolutionary 'shift' events that allow some incongruence between sequence and structure evolution.

![](Examples/example.svg)

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

## Usage examples

The tool can be used from the command line or via its Python interface
(e.g. from a Jupyter notebook).
### Command line interface

To get an overview on all command line parameters that configure modes and
alignment parameters, please refer to the help output of the tool as obtained by

```
./bialign.py --help
```

#### RNA bi-alignment examples

This 'toy' example demonstrates a simple helix shift;

```bash
./bialign.py GCGGGGGAUAUCCCCAUCG GGGGAUAUCCCCAUCG \
    --strA "...(((.....)))....." --strB ".(((.....)))...." \
    --structure 400 \
    --gap_opening_cost -200 --gap_cost -50 \
    --max_shift 0 --shift_cost -150
```

Structures will be predicted (using the Vienna RNA package) if they are not
explicitly given, e.g.
```
./bialign.py UGUAAACAUCCUCGACUGGAAGCUGUGAAGCCACAAAUGGGCUUUCAGUCGGAUGUUUGCA UGUAAACAUCCUACACUCAGCUGUCAUACAUGCGUUGGCUGGGAUGUGGAUGUUUACG 
```
Note that this fails, if the Vienna RNA package with Python binding is not
available.


#### Bi-Alignments of proteins with affine gap cost

```bash
./bialign.py RAKLPLKEKKLTATANYHPGIRYIMTGYSAKYIYSSTYARFR KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAYFR \
    --strA "CHHHHHHHHHHHHHCCCCTCEEEEEEECCTCEEEEEEEECCC" --strB "HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEECC" \
    --type Protein --shift_cost -150 --structure_weight 800 --simmatrix BLOSUM62 --gap_opening_cost -150 \
    --gap_cost -50 --max_shift 1
```

Input can also be read from files as written by the secondary structure prediction
web server CFSSP (Kumar et al, 2013; http://www.biogem.org/tool/chou-fasman).

```bash
./bialign.py --filein Examples/DNAPolymerase1_Escherichia.cfssp Examples/DNAPolymerase1_Xanthomonas.cfssp \
    --type Protein --shift_cost -150 --structure_weight 800 --simmatrix BLOSUM62 --gap_opening_cost -150 \
    --gap_cost -50 --max_shift 1
```


### Python interface

The following code generates a bi-alignment of two toy proteins and shows
the resulting alignment in a graphical representation.

```python
import bialignment
import bialignment as ba
import timeit

args = {'type': 'Protein',
        'gap_cost': -50,
        'gap_opening_cost': -150,
        'shift_cost': -150,
        'structure_weight': 800,
        'max_shift': 1,
        'simmatrix': 'BLOSUM62'
       }

args['nameA'] = 'A'
args['nameB'] = 'B'
strA = "CHHHHHHHHHHHHHCCCCTCEEEEEEECCTCEEEEEEEECCC"
seqA = "RAKLPLKEKKLTATANYHPGIRYIMTGYSAKYIYSSTYARFR"
seqB = "KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAYFR"
strB = "HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEECC"


bialigner = ba.BiAligner(seqA, seqB, strA, strB,
                         **args)

score = bialigner.optimize()
print('SCORE',score)
print()

for line in bialigner.decode_trace():
    print(line)

ba.plot_alignment(bialigner.decode_trace_full(),
    width = 80,
    show_position_numbers=False,
    name_offset=3,
    #outname = "example.svg" #optionally write plot to file
    )
```
![](Examples/example.svg)
