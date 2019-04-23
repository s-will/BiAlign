#!/usr/bin/env python3
import RNA

import sys
import argparse
import numpy as np
from math import log,exp,sqrt

## Alignment factory
class BiAligner:
    def __init__(self, seqA, seqB, strA=None, strB=None):
        self.rnaA = self._preprocess_seq(seqA,strA)
        self.rnaB = self._preprocess_seq(seqB,strB)

        # parametrization
        self._sequence_match_similarity = 102
        self._sequence_mismatch_similarity = -50
        self._structure_weight = 100
        self._gap_cost = -50

        self._shift_cost = -50 # cost of shifting the 2 scores against each other
        self._max_shift = 3 # maximal number of shifts away from the diagonal in either direction

        print( "Structure similarity",
               [ ((i,j),s) for i in range(1,len(seqA)+1)
                           for j in range(1,len(seqB)+1)
                           for s in [self._structure_similarity(i,j)] if s!=0])

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
        yield ((1,1,0,0), self.mu1(i,j) + self._shift_cost)
        yield ((0,0,1,1), self.mu2(i,j) + self._shift_cost)

        yield ((1,0,0,0), self.g1A(i) + self._shift_cost)
        yield ((0,1,0,0), self.g1B(j) + self._shift_cost)
        yield ((0,0,1,0), self.g2A(k) + self._shift_cost)
        yield ((0,0,0,1), self.g2B(l) + self._shift_cost)

        yield ((1,0,1,1), self.g1A(i) + self.mu2(k,l) + self._shift_cost)
        yield ((0,1,1,1), self.g1B(j) + self.mu2(k,l) + self._shift_cost)
        yield ((1,1,1,0), self.g2A(k) + self.mu1(i,j) + self._shift_cost)
        yield ((1,1,0,1), self.g2B(l) + self.mu1(i,j) + self._shift_cost)

        # double-shift cases -- these cases can be replaced by two others --> skip
        # yield ((0,1,1,0), self.g1B(j) + self.g2A(k) + 2 * self._shift_cost)
        # yield ((1,0,0,1), self.g1A(i) + self.g2B(l) + 2 * self._shift_cost)

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
        return i-io>=0 and j-jo>=0 and k-ko>=0 and l-lo >=0 and abs(k-ko-(i-io))<=self._max_shift and abs(l-lo-(j-jo))<=self._max_shift

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
        else:
            if len(structure)!=len(sequence):
                print("Fixed structure and sequence must have the same length.")
                sys.exit()
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
            return self._sequence_match_similarity
        else:
            return self._sequence_mismatch_similarity

    def _structure_similarity(self,i,j):
        return int( self._structure_weight *
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
        return self._gap_cost
    # gap cost for inserting i, score 1, in second sequence
    def g1B(self,i):
        return self._gap_cost
    # gap cost for inserting i, score 2, in first sequence
    def g2A(self,i):
        return self._gap_cost
    # gap cost for inserting i, score 2, in second sequence
    def g2B(self,i):
        return self._gap_cost

    # run alignment algorithm
    def optimize(self):
        lenA = self.rnaA["len"]
        lenB = self.rnaB["len"]

        self.M = np.zeros((lenA+1,lenB+1,lenA+1,lenB+1), dtype=int)

        for i in range(0,lenA+1):
            for j in range(0,lenB+1):
                for k in range( max(0, i-self._max_shift), min(lenA+1, i+self._max_shift+1) ):
                    for l in range( max(0, j-self._max_shift), min(lenB+1, j+self._max_shift+1) ):
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

    # decode trace to alignment strings
    def decode_trace(self,trace):
        seqs = (self.rnaA["seq"],self.rnaB["seq"],self.rnaA["seq"],self.rnaB["seq"])
        pos = [0]*4
        alignment = [""]*4
        for i,y in enumerate(trace):
            for s in range(4):
                if (y[s]==0):
                    alignment[s] = alignment[s] + "-"
                elif (y[s]==1):
                    alignment[s] = alignment[s] + seqs[s][pos[s]]
                    pos[s]+=1
        return alignment

    # decode trace to alignment strings
    def eval_trace(self,trace):
        pos=[0]*4
        for i,y in enumerate(trace):
            for k in range(4): pos[k] += y[k]
            # lookup case
            for x in self.recursionCases(pos[0],pos[1],pos[2],pos[3]):
                if x[0] == y:
                    print(pos,
                          y,
                          x[1],
                          "-->",
                          self.evalCase(x,pos[0],pos[1],pos[2],pos[3]))
                    break


def main(args):
    ba = BiAligner(args.seqA,args.seqB,args.strA,args.strB)

    optscore = ba.optimize()
    print("SCORE:",optscore)
    trace    = ba.traceback()
    for s in ba.decode_trace(trace): print(s)
    if args.verbose:
        ba.eval_trace(trace)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= "Bialignment.")
    parser.add_argument("seqA",help="RNA sequence A")
    parser.add_argument("seqB",help="RNA sequence B")
    parser.add_argument("--strA",default=None,help="RNA structure A")
    parser.add_argument("--strB",default=None,help="RNA structure B")
    parser.add_argument("-v","--verbose",action='store_true',help="Verbose")

    args = parser.parse_args()
    main(args)
