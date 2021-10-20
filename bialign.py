#!/usr/bin/env python3
import argparse

import bialignment as ba
import sys
from collections import defaultdict

VERSION_STRING = f"BiAlign {ba.__version__}"

def read_molecule(content):
    result = defaultdict(lambda:"")
    keys =["Query","Struc"]
    for line in content.split('\n'):
        line = line.split()
        if not line:
            continue
        if line[0] in keys:
            if len(line) != 4:
                raise IOError("Cannot parse")
            result[line[0]] += line[2]

    if len(result[keys[0]]) != len(result[keys[1]]):
        raise IOError("Sequence and structure of unequal length.")
    if len(result[keys[0]]) == 0:
        raise IOError("Input does not contain input sequence and structure.")

    return [result[k] for k in keys]

def read_molecule_from_file(filename):
    try:
        with open(filename,'r') as fh:
            return read_molecule(fh.read())
    except FileNotFoundError as e:
        print("Input file not found.")
        print(e)
        sys.exit(-1)
    except IOError as e:
        print(f"Cannot read input file {filename}.")
        print(e)
        sys.exit(-1)

def add_bialign_parameters(parser):
    parser.add_argument("seqA", help="sequence A")
    parser.add_argument("seqB", help="sequence B")
    parser.add_argument("--strA", default=None, help="structure A")
    parser.add_argument("--strB", default=None, help="structure B")
    parser.add_argument("--nameA", default="A", help="name A")
    parser.add_argument("--nameB", default="B", help="name B")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose")

    parser.add_argument(
        "--type", default="RNA", type=str, help="Type of molecule: RNA or Protein"
    )

    parser.add_argument(
        "--nodescription",
        action="store_true",
        help="Don't prefix the strings in output alignment with descriptions",
    )
    parser.add_argument(
        "--outmode",
        default="sorted",
        help="Output mode [call --mode help for a list of options]",
    )

    parser.add_argument(
        "--sequence_match_similarity",
        type=int,
        default=100,
        help="Similarity of matching nucleotides",
    )
    parser.add_argument(
        "--sequence_mismatch_similarity",
        type=int,
        default=0,
        help="Similarity of mismatching nucleotides",
    )
    parser.add_argument(
        "--structure_weight",
        type=int,
        default=100,
        help="Weighting factor for structure similarity",
    )
    parser.add_argument(
        "--gap_opening_cost",
        type=int,
        default=0,
        help="Similarity of opening a gap (turns on affine gap cost if not 0)",
    )
    parser.add_argument(
        "--gap_cost", type=int, default=-200, help="Similarity of a single gap position"
    )
    parser.add_argument(
        "--shift_cost",
        type=int,
        default=-250,
        help="Similarity of shifting the two scores against each other",
    )
    parser.add_argument(
        "--max_shift",
        type=int,
        default=2,
        help="Maximal number of shifts away from the diagonal in either direction",
    )
    parser.add_argument(
        "--fileinput",
        action="store_true",
        help="Read sequence and structure input from file",
    )

    parser.add_argument("--version", action="version", version=VERSION_STRING)

    parser.add_argument("--simmatrix", type=str, default=None, help="Similarity matrix")


def main():
    parser = argparse.ArgumentParser(description="Bialignment.")
    add_bialign_parameters(parser)

    args = parser.parse_args()

    if args.seqA == "example":
        args.seqA = "MSKLVLIDGSSYLYRAFHALPPLTNAQGEPTGALFGVVNMLRATLKERPAYVAFVVDAPGKTFRDDLYADYKANRPSMPDELRAQVQPMCDIVHALGIDILRIDGVEADD"
        args.strA = "HEEEEEHCTTCEEEEHHCCCCCCCCCCTTCCCHEEEEEHHHHHHHHHTTHEEEEEHHCCTTCCCTCCCCCCCCCTTCCCHHHHEEEEHEEEEEHEEEEEEEHHHHHHHHH"
        args.seqB = "MVQIPQNPLILVDGSSYLYRAYHAFPPLTNSAGEPTGAMYGVLNMLRSLIMQYKPTHAAVVFDAKGKTFRDELFEHYKSHRPPMPDDLRAQIEPLHAMVKAMGLPLLAVS"
        args.strB = "EEEEETEEEEEHCTTCEEEEEECCCCCCCTCCCTCCCHEEEEEHHEEEEEHEHCTCHHHHHHHHHTHHHHHHHHHHHHTCCTCCCTHHHHHHHHHHHHHHHHEEHEEEEH"

        L = 60
        args.seqA = args.seqA[:L]
        args.strA = args.strA[:L]
        args.seqB = args.seqB[:L]
        args.strB = args.strB[:L]

    if args.fileinput:
        args.seqA, args.strA = read_molecule_from_file(args.seqA)
        args.seqB, args.strB = read_molecule_from_file(args.seqB)


    print(
        "\n".join(
            [
                "Input:",
                "seqA\t " + args.seqA,
                "seqB\t " + args.seqB,
                "strA\t " + args.strA,
                "strB\t " + args.strB,
                "",
            ]
        )
    )

    if args.outmode == "help":
        print()
        print("Available modes: " + ", ".join(BiAligner.outmodes.keys()))
        print()
        exit()

    for line in ba.bialign(**vars(args)):
        print(line)


if __name__ == "__main__":
    main()

