# BiAlign - Bialignment of RNAs and proteins

The tool BiAlign computes optimal bi-alignments of RNAs and proteins. Such
bi-alignments support evolutionary 'shift' events between sequence and
structure. In this way, bialignments extend alignments based on sequence
and struture similarity to the case of potential incongruence between
sequence and structure evolution.

The current version extends the capabilities from RNA alignments to the alignment of protein
sequence and secondary structure, supporting realistic 'affine' gap cost
with gap opening and extension scores.

![](Examples/example.svg)

The first version of this tool has been described in

Waldl M., Will S., Wolfinger M.T., Hofacker I.L., Stadler P.F. (2020)
Bi-alignments as Models of Incongruent Evolution of RNA Sequence and
Secondary Structure. In: Cazzaniga P., Besozzi D., Merelli I., Manzoni L.
(eds) Computational Intelligence Methods for Bioinformatics and
Biostatistics. CIBB 2019. Lecture Notes in Computer Science, vol 12313.
Springer, Cham. https://doi.org/10.1007/978-3-030-63061-4_15



## Installation

This software will run with full functionality only on Linux and Mac
systems. Installation via conda is not supported on Windows
and the prediction of RNA structures (using the Vienna RNA package) cannot be
supported.

The software can be installed via Conda (only Linux/Mac) or pip (Mac/Linux/Windows/...)  respectively by

```
conda install -c bioconda bialign
```

or

```
pip install bialign
```


Conda installation is recommended, since it will automatically install dependencies like the Vienna RNA
package. When installing via pip (or from source, see below), additionally install numpy, matplotlib, and (optionally) the Vienna RNA package.


### Installation from source
Installation or from source, e.g. a clone of the git repository, relies on
the python setup system.

We require Cython to compile performance critical code. For this purpose Cython and
Python (including pip/setuptools) need to be installed. Install from source by

```
pip install .
```

Moreover, for aligning RNAs, the tool requires the Vienna RNA package with Python bindings. We recommend to use the tool under Linux or MacOS and install the prerequisites via conda / bioconda.

## Usage examples

The tool can be used from the command line or via its Python interface
(e.g. from a Jupyter notebook).

### Command line interface

To get an overview on all command line parameters that configure modes and
alignment parameters, please refer to the help output of the tool as obtained by

```
bialign.py --help
```

#### RNA bi-alignment examples

This 'toy' example demonstrates a simple helix shift:

```bash
bialign.py GCGGGGGAUAUCCCCAUCG GGGGAUAUCCCCAUCG \
    --strA "...(((.....)))....." --strB ".(((.....)))...." \
    --structure 400 \
    --gap_opening_cost -200 --gap_cost -50 \
    --max_shift 1 --shift_cost -150
```

Using default text output mode, this produces
```
Input:
seqA	 GCGGGGGAUAUCCCCAUCG
seqB	 GGGGAUAUCCCCAUCG
strA	 ...(((.....))).....
strB	 .(((.....)))....
SCORE: 6800

A               GCGGGGGAUAUCCCC-AUCG
B               G---GGGAUAUCCCC-AUCG
A ss            ...-(((.....))).....
B ss            .---(((.....)))-....
A shifts        ...<...........>....
B shifts        ....................
```

Structures will be predicted (using the Vienna RNA package) if they are not
explicitly given, e.g.
```
bialign.py UGUAAACAUCCUCGACUGGAAGCUGUGAAGCCACAAAUGGGCUUUCAGUCGGAUGUUUGCA UGUAAACAUCCUACACUCAGCUGUCAUACAUGCGUUGGCUGGGAUGUGGAUGUUUACG
```
Note that this fails, if the Vienna RNA package with Python binding is not
available.


#### Bi-Alignments of proteins with affine gap cost

```bash
bialign.py RAKLPLKEKKLTATANYHPGIRYIMTGYSAKYIYSSTYARFR KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAYFR \
    --strA "CHHHHHHHHHHHHHCCCCTCEEEEEEECCTCEEEEEEEECCC" --strB "HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEECC" \
    --type Protein --shift_cost -150 --structure_weight 800 --simmatrix BLOSUM62 --gap_opening_cost -150 \
    --gap_cost -50 --max_shift 1 --outmode sorted
```

Due to the requested output mode `sorted`, this produces text output with BLAST-like
annotation by the respective consensus sequence and structure of
the sequence and structure alignment component.

```
Input:
seqA	 RAKLPLKEKKLTATANYHPGIRYIMTGYSAKYIYSSTYARFR
seqB	 KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAYFR
strA	 CHHHHHHHHHHHHHCCCCTCEEEEEEECCTCEEEEEEEECCC
strB	 HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEECC
SCORE: 48500

A ss            -CHHHHHHHHHHHHHCCCCTCEEEEEEECCTCEEEEEEEEC-CC
A               -RAKLPLKEKKLTATANYHPGIRYIMTGYSAKYIYSSTYAR-FR
consensus       -.AKLPLKEKKLT.TANYHPGIRYIMTGYSAK.IYSSTYA.-FR
B               -KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAY-FR
B ss            -HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEE-CC
consensus ss    -.HHHHHHHHHHH..CCCCTCEEEEEEECC.C.EEEEEEE.-CC

A               RAKLPLKEKKLTA-TANYHPGIRYIMTGYSAK-YIYSSTYARFR
A ss            CHHHHHHHHHHHH-HCCCCTCEEEEEEECCTC-EEEEEEEECCC
consensus ss    .HHHHHHHHHHHH..CCCCTCEEEEEEECC.C.EEEEEEEE.CC
B ss            -HHHHHHHHHHHHCCCCCCTCEEEEEEECCCCCEEEEEEEE-CC
B               -KAKLPLKEKKLTRTANYHPGIRYIMTGYSAKRIYSSTYAY-FR
consensus       .........K....TANYHPGIRYIMTGYSAK....S.....FR

A shifts        >............<..................<........>..
B shifts        ............................................
```


Input can also be read from files as written by the secondary structure prediction
web server CFSSP (Kumar et al, 2013; http://www.biogem.org/tool/chou-fasman).

```bash
bialign.py --filein Examples/DNAPolymerase1_Escherichia.cfssp Examples/DNAPolymerase1_Xanthomonas.cfssp \
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
