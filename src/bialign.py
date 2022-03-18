#!/usr/bin/env python3
import argparse

import bialignment
import sys

VERSION_STRING = f"BiAlign {bialignment.__version__}"


def bialign(seqA, seqB, strA, strB, verbose, **args):
    ba = bialignment.BiAligner(seqA, seqB, strA, strB, **args)

    optscore = ba.optimize()
    yield "SCORE: " + str(optscore)
    yield ""

    ali = ba.decode_trace()

    yield from ali

    if verbose:
        yield from ba.eval_trace()


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
        default="default",
        help="Output mode [call --outmode help for a list of options]",
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
        default=400,
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

    if args.fileinput:
        args.seqA, args.strA = bialignment.read_molecule_from_file(args.seqA, args.type)
        args.seqB, args.strB = bialignment.read_molecule_from_file(args.seqB, args.type)

    input_descr = ["Input:", "seqA\t " + args.seqA, "seqB\t " + args.seqB]
    if hasattr(args, "strA") and args.strA is not None:
        input_descr.append("strA\t " + args.strA)
    if hasattr(args, "strB") and args.strB is not None:
        input_descr.append("strB\t " + args.strB)

    print("\n".join(input_descr))

    if args.outmode == "help":
        print()
        print("Available modes: " + ", ".join(bialignment.BiAligner.outmodes.keys()))
        print()
        exit()

    for line in bialign(**vars(args)):
        print(line)


if __name__ == "__main__":
    main()
