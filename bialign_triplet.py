#!/usr/bin/env python3
import RNA

import sys
import argparse
import numpy as np
from math import log,exp

from bialign import BiAligner

## Alignment factory
class BiAlignerTriplet (BiAligner):
    def __init__(self, seqA, seqB, strA=None, strB=None):
        super().__init__(seqA, seqB, strA, strB)

    # iterator over recursion cases
    #
    # per case:
    #       pair of access info
    #     and
    #       function that returns list of case score components
    def recursionCases(self,i,j,k):
        # synchronous cases
        yield ((1,1,1), self.mu1(i,j) + self.mu2(i,k))
        yield ((1,0,0), self.g1A(i)   + self.g2A(i))
        yield ((0,1,1), self.g1B(j)   + self.g2B(k))
        # shifting
        yield ((1,1,0), self.mu1(i,j) + self.g2A(i) + self._shift_cost)
        yield ((1,0,1), self.mu2(i,k) + self.g1A(i) + self._shift_cost)
        yield ((0,1,0), self.g1B(j)   + self._shift_cost)
        yield ((0,0,1), self.g2B(k)   + self._shift_cost)

    def guardCase(self,x,i,j,k):
        (io,jo,ko) = x[0]
        return i-io>=0 and j-jo>=0 and k-ko>=0 and abs(k-ko-(j-jo))<=self._max_shift

    def evalCase(self,x,i,j,k):
        (io,jo,ko) = x[0]
        return self.M[i-io,j-jo,k-ko] + x[1]

    # run alignment algorithm
    def optimize(self):
        lenA = self.rnaA["len"]
        lenB = self.rnaB["len"]

        self.M = np.zeros((lenA+1,lenB+1,lenB+1), dtype=int)

        for i in range(0,lenA+1):
            for j in range(0,lenB+1):
                for k in range( max(0,j-self._max_shift), min(lenB+1, j+self._max_shift+1) ):
                    self.M[i,j,k] = self.plus( [ self.evalCase(x,i,j,k)
                                                 for x in self.recursionCases(i,j,k)
                                                 if self.guardCase(x,i,j,k) ] )
        return self.M[lenA,lenB,lenB]

    # perform traceback
    # @returns list of 'trace arrows'
    def traceback(self):
        lenA = self.rnaA["len"]
        lenB = self.rnaB["len"]

        trace=[]

        def trace_from(i,j,k):
            for x in self.recursionCases(i,j,k):
                if self.guardCase(x,i,j,k):
                    if self.evalCase(x,i,j,k) == self.M[i,j,k]:
                        (io,jo,ko) = x[0]
                        trace.append((io,jo,ko))
                        trace_from(i-io,j-jo,k-ko)
                        break

        trace_from(lenA,lenB,lenB)
        return list(reversed(trace))

    # decode trace to alignment strings
    def decode_trace(self,trace):
        seqs = (self.rnaA["seq"],self.rnaB["seq"],self.rnaB["seq"])
        pos = [0]*3
        alignment = [""]*3
        for i,y in enumerate(trace):
            for s in range(3):
                if (y[s]==0):
                    alignment[s] = alignment[s] + "-"
                elif (y[s]==1):
                    alignment[s] = alignment[s] + seqs[s][pos[s]]
                    pos[s]+=1
        return alignment

    # decode trace to alignment strings
    def eval_trace(self,trace):
        pos=[0]*3
        for i,y in enumerate(trace):
            for k in range(3): pos[k] += y[k]
            # lookup case
            for x in self.recursionCases(pos[0],pos[1],pos[2]):
                if x[0] == y:
                    print(pos,
                          y,
                          x[1],
                          "-->",
                          self.evalCase(x,pos[0],pos[1],pos[2]))
                    break


def main(args):
    ba = BiAlignerTriplet(args.seqA,args.seqB,args.strA,args.strB)

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
