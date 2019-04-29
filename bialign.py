#!/usr/bin/env python3

import RNA

import sys
import argparse
import numpy as np
from math import log,exp,sqrt

## Alignment factory
class BiAligner:
    def __init__(self, seqA, seqB, strA, strB, **params):
        self.rnaA = self._preprocess_seq(seqA,strA)
        self.rnaB = self._preprocess_seq(seqB,strB)

        # parametrization
        self._params = params

        #print( "Structure similarity",
        #       [ ((i,j),s) for i in range(1,len(seqA)+1)
        #                   for j in range(1,len(seqB)+1)
        #                   for s in [self._structure_similarity(i,j)] if s!=0])

        # the dynamic programming matrix
        self.M = None

    # GENERAL ALGORITHM IDEAS
    #
    # - perform bialignment
    # - use two copies of the input sequences
    # - use indices i and j for the first copy; k and l for the second
    # - score and restrict shift


    # iterator over recursion cases
    #
    # per case:
    #       pair of access info
    #     and
    #       function that returns list of case score components
    def recursionCases(self,i,j,k,l):
        # synchronous cases
        yield ((1,1,1,1), self.mu1(i,j) + self.mu2(k,l))
        yield ((1,0,1,0), self.g1A(i)   + self.g2A(k))
        yield ((0,1,0,1), self.g1B(j)   + self.g2B(l))
        # shifting
        yield ((1,1,0,0), self.mu1(i,j) + self._params["shift_cost"])
        yield ((0,0,1,1), self.mu2(i,j) + self._params["shift_cost"])

        yield ((1,0,0,0), self.g1A(i) + self._params["shift_cost"])
        yield ((0,1,0,0), self.g1B(j) + self._params["shift_cost"])
        yield ((0,0,1,0), self.g2A(k) + self._params["shift_cost"])
        yield ((0,0,0,1), self.g2B(l) + self._params["shift_cost"])

        yield ((1,0,1,1), self.g1A(i) + self.mu2(k,l) + self._params["shift_cost"])
        yield ((0,1,1,1), self.g1B(j) + self.mu2(k,l) + self._params["shift_cost"])
        yield ((1,1,1,0), self.g2A(k) + self.mu1(i,j) + self._params["shift_cost"])
        yield ((1,1,0,1), self.g2B(l) + self.mu1(i,j) + self._params["shift_cost"])

        # double-shift cases -- these cases can be replaced by two others --> skip
        # yield ((0,1,1,0), self.g1B(j) + self.g2A(k) + 2 * self._params["shift_cost"])
        # yield ((1,0,0,1), self.g1A(i) + self.g2B(l) + 2 * self._params["shift_cost"])

    # plus operator (max in optimization; sum in pf)
    def plus(self, xs):
        if xs != []:
           return max(xs)
        else:
           return 0

    ## mul operator (sum in optimization, product in pf)
    # def mul(self, xs):
    #    return sum(xs)

    def guardCase(self,x,i,j,k,l):
        (io,jo,ko,lo) = x[0]
        return i-io>=0 and j-jo>=0 and k-ko>=0 and l-lo >=0 and abs(k-ko-(i-io))<=self._params["max_shift"] and abs(l-lo-(j-jo))<=self._params["max_shift"]

    def evalCase(self,x,i,j,k,l):
        (io,jo,ko,lo) = x[0]
        return self.M[i-io,j-jo,k-ko,l-lo] + x[1]

    # make bpp symmetric (based on upper triangular matrix)
    # and set diagonal to unpaired probs
    @staticmethod
    def _symmetrize_bpps(bpp):
        n=len(bpp)-1
        sbpp = np.zeros((n+1,n+1),dtype='float')
        for i in range(1,n+1):
            for j in range(i+1,n+1):
                sbpp[i,j] = bpp[i][j]
                sbpp[j,i] = bpp[i][j]

        for i in range(1,n+1):
            sbpp[i,i] = 1.0 - sum( sbpp[i,j] for j in range(1,n+1) )

        return sbpp

    @staticmethod
    def _preprocess_seq(sequence, structure):
        x = dict()
        x["seq"] = str(sequence)
        x["len"] = len(x["seq"])

        if structure is None:
            fc = RNA.fold_compound(str(sequence))
            x["mfe"] = fc.mfe()
            x["pf"] = fc.pf()
            x["sbpp"] = BiAligner._symmetrize_bpps( fc.bpp() )
            x["structure"] = x["mfe"][0]
        else:
            if len(structure)!=len(sequence):
                print("Fixed structure and sequence must have the same length.")
                sys.exit()
            x["structure"] = structure
            x["sbpp"] = BiAligner._bp_matrix_from_fixed_structure(structure)

        n = x["len"]
        #note: the arrays are 1-based (we just ignore entries at 0)
        x["up"]   = [ sum( x["sbpp"][i][j] for j in range(1,i-1) ) for i in range(0,n+1) ]
        x["down"] = [ sum( x["sbpp"][i][j] for j in range(i+1,n+1) ) for i in range(0,n+1) ]
        x["unp"]  = [ 1.0 - x["up"][i] - x["down"][i] for i in range(0,n+1) ]

        return x

    @staticmethod
    def _bp_matrix_from_fixed_structure(structure):
        n = len(structure)
        bpm = np.zeros((n+1,n+1),dtype='float')
        stack=list()
        for i in range(n):
            if structure[i]=='(':
                stack.append(i)
            elif structure[i]==')':
                j=stack.pop()
                bpm[i+1,j+1]=1.0
                bpm[j+1,i+1]=1.0
            else:
                bpm[i+1,i+1]=1.0
        return bpm

    @staticmethod
    def _expected_pairing(rna):
        n = rna["len"]
        sbpp = rna["sbpp"]
        def ep(i):
            return sum( sbpp[i,j]*(j-i) for j in range(1,n+1) )

        return [0] + [ ep(i) for i in range(1,n+1) ]

    # sequence similarity of residues i and j, 1-based
    def _sequence_similarity(self,i,j):
        if self.rnaA["seq"][i-1]==self.rnaB["seq"][j-1]:
            return self._params["sequence_match_similarity"]
        else:
            return self._params["sequence_mismatch_similarity"]

    def _structure_similarity(self,i,j):
        return int( self._params["structure_weight"] *
                    (
                      sqrt(self.rnaA["up"][i]*self.rnaB["up"][j])
                      + sqrt(self.rnaA["down"][i]*self.rnaB["down"][j])
                      + sqrt(self.rnaA["unp"][i]*self.rnaB["unp"][j])
                    )
                  )
    # Scoring functions
    # note: scoring functions have 1-based indices

    # match/mismatch cost for i~j, score 1
    def mu1(self,i,j):
        return self._sequence_similarity(i,j)

    # match/mismatch cost for i~j, score 2
    def mu2(self,i,j):
        return self._structure_similarity(i,j)

    # gap cost for inserting i, score 1, in first sequence
    def g1A(self,i):
        return self._params["gap_cost"]
    # gap cost for inserting i, score 1, in second sequence
    def g1B(self,i):
        return self._params["gap_cost"]
    # gap cost for inserting i, score 2, in first sequence
    def g2A(self,i):
        return self._params["gap_cost"]
    # gap cost for inserting i, score 2, in second sequence
    def g2B(self,i):
        return self._params["gap_cost"]

    # run alignment algorithm
    def optimize(self):
        lenA = self.rnaA["len"]
        lenB = self.rnaB["len"]

        self.M = np.zeros((lenA+1,lenB+1,lenA+1,lenB+1), dtype=int)

        for i in range(0,lenA+1):
            for j in range(0,lenB+1):
                for k in range( max(0, i-self._params["max_shift"]), min(lenA+1, i+self._params["max_shift"]+1) ):
                    for l in range( max(0, j-self._params["max_shift"]), min(lenB+1, j+self._params["max_shift"]+1) ):
                        self.M[i,j,k,l] = self.plus( [ self.evalCase(x,i,j,k,l)
                                                     for x in self.recursionCases(i,j,k,l)
                                                     if self.guardCase(x,i,j,k,l) ] )
        return self.M[lenA,lenB,lenA,lenB]

    # perform traceback
    # @returns list of 'trace arrows'
    def traceback(self):
        lenA = self.rnaA["len"]
        lenB = self.rnaB["len"]

        trace=[]

        def trace_from(i,j,k,l):
            for x in self.recursionCases(i,j,k,l):
                if self.guardCase(x,i,j,k,l):
                    if self.evalCase(x,i,j,k,l) == self.M[i,j,k,l]:
                        (io,jo,ko,lo) = x[0]
                        trace.append((io,jo,ko,lo))
                        trace_from(i-io,j-jo,k-ko,l-lo)
                        break

        trace_from(lenA,lenB,lenA,lenB)
        return list(reversed(trace))

    # transfer gap pattern from an alignment string to a sequence string
    @staticmethod
    def _transfer_gaps(alistr,seqstr):
        pos=0
        res=""
        for i,c in enumerate(alistr):
            if c=="-":
                res+="-"
            else:
                res+=seqstr[pos]
                pos+=1
        return res

    # decode trace to alignment strings
    def decode_trace(self, *, trace=None,
            show_structures=False,highlight_identity=False):
        if trace is None:
            trace = self.traceback()

        rnas = (self.rnaA,self.rnaB,self.rnaA,self.rnaB)
        pos = [0]*len(rnas)
        alignment = [""]*len(rnas)
        for i,y in enumerate(trace):
            for s in range(len(rnas)):
                if (y[s]==0):
                    alignment[s] = alignment[s] + "-"
                elif (y[s]==1):
                    alignment[s] = alignment[s] + rnas[s]["seq"][pos[s]]
                    pos[s]+=1

        if highlight_identity:
            alignment[0],alignment[1] = highlight_sequence_identity(alignment[0],alignment[1])
            alignment[2],alignment[3] = highlight_sequence_identity(alignment[2],alignment[3])

        if not show_structures:
            return alignment


        # annotate with structure
        anno_ali = list()
        for alistr,rna in zip(alignment,rnas):
            anno_ali.append( self._transfer_gaps(alistr,rna["structure"]) )
            anno_ali.append( alistr )

        if highlight_identity:
            anno_ali[0],anno_ali[2] = highlight_structure_identity(anno_ali[0],anno_ali[2])
            anno_ali[4],anno_ali[6] = highlight_structure_identity(anno_ali[4],anno_ali[6])

        return anno_ali

    # evaluate trace
    def eval_trace(self, trace=None):
        if trace is None:
            trace = self.traceback()

        pos=[0]*4
        for i,y in enumerate(trace):
            for k in range(4): pos[k] += y[k]
            # lookup case
            for x in self.recursionCases(pos[0],pos[1],pos[2],pos[3]):
                if x[0] == y:
                    line = " ".join([ str(x) for x in
                            [pos,
                            y,
                            x[1],
                            "-->",
                            self.evalCase(x,pos[0],pos[1],pos[2],pos[3])
                            ]
                        ])
                    yield line
                    break

# highlight identical sequence in a pairwise alignment
def highlight_sequence_identity(alistrA,alistrB):
    res = ["",""]
    for x,y in zip(alistrA.lower(),alistrB.lower()):
        if x==y:
            x=x.upper()
            y=x
        res[0]+=x
        res[1]+=y
    return res


def parse_dotbracket(dbstr):
    res = [-1] * len(dbstr)
    stack=list()
    for i,sym in enumerate(dbstr):
        if sym=='(':
            stack.append(i)
        elif sym==')':
            j=stack.pop()
            res[i]=j
            res[j]=i

    return res

# highlight matched base pairs in a pairwise alignment
def highlight_structure_identity(alistrA,alistrB):

    strA = parse_dotbracket(alistrA)
    strB = parse_dotbracket(alistrB)

    res = ["",""]
    for i,(x,y) in enumerate(zip(alistrA.lower(),alistrB.lower())):
        if strA[i]>=0 and strB[i]>=0 and strA[i]==strB[i]:
            if strA[i]>i:
                x='['
            else:
                x=']'
            y=x

        res[0]+=x
        res[1]+=y

    return res

def bialign(seqA,seqB,strA,strB,show_structures,highlight_identity,verbose,**args):
    ba = BiAligner(seqA,seqB,strA,strB,**args)

    optscore = ba.optimize()
    yield "SCORE:" + str(optscore)

    ali = ba.decode_trace(show_structures=show_structures,highlight_identity=highlight_identity)

    yield from ali

    if verbose:
        yield from ba.eval_trace()

def add_bialign_parameters(parser):
    parser.add_argument("seqA",help="RNA sequence A")
    parser.add_argument("seqB",help="RNA sequence B")
    parser.add_argument("--strA",default=None,help="RNA structure A")
    parser.add_argument("--strB",default=None,help="RNA structure B")
    parser.add_argument("--show_structures",action='store_true',help="Print structure annotation with alignment")
    parser.add_argument("-v","--verbose",action='store_true',help="Verbose")

    parser.add_argument("--highlight_identity",action='store_true',help="Highlight identity")

    parser.add_argument("--sequence_match_similarity", type=int,
            default=100, help="Similarity of matching nucleotides")
    parser.add_argument("--sequence_mismatch_similarity", type=int,
            default=-50, help="Similarity of mismatching nucleotides")
    parser.add_argument("--structure_weight", type=int, default=100,
            help="Weighting factor for structure similarity")
    parser.add_argument("--gap_cost", type=int, default=-50,
            help="Similarity of a single gap")
    parser.add_argument("--shift_cost", type=int, default=-50,
            help="Similarity of shifting the two scores against each other")
    parser.add_argument("--max_shift", type=int, default=3,
            help="Maximal number of shifts away from the diagonal in either direction")

def main():
    parser = argparse.ArgumentParser(description= "Bialignment.")
    add_bialign_parameters(parser)

    args = parser.parse_args()

    for line in bialign(**vars(args)):
        print(line)

if __name__ == "__main__":
    main()
