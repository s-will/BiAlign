#!/usr/bin/env python3
import RNA

import sys
import argparse
import numpy as np
from math import log,exp

from bialign import *

## Alignment factory
class BiAlignerTriplet (BiAligner):
    def __init__(self, seqA, seqB, strA, strB, **params):
        super().__init__(seqA, seqB, strA, strB, **params)

    # iterator over recursion cases
    #
    # per case:
    #       pair of access info
    #     and
    #       function that returns list of case score components
    def recursionCases(self,i,j,k):
        mu1ij = self.mu1(i,j)
        mu2ik = self.mu2(i,k)
        Delta = self._params["shift_cost"]

        # synchronous cases
        yield ((1,1,1), mu1ij + mu2ik))
        yield ((1,0,0), self.g1A(i) + self.g2A(i))
        yield ((0,1,1), self.g1B(j) + self.g2B(k))
        # shifting
        yield ((1,1,0), mu1ij + self.g2A(i) + Delta)
        yield ((1,0,1), mu2ik + self.g1A(i) + Delta)
        yield ((0,1,0), self.g1B(j) + Delta)
        yield ((0,0,1), self.g2B(k) + Delta)

    def guardCase(self,x,i,j,k):
        (io,jo,ko) = x[0]
        return i-io>=0 and j-jo>=0 and k-ko>=0 and abs(k-ko-(j-jo))<=self._params["max_shift"]

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
                for k in range( max(0,j-self._params["max_shift"]), min(lenB+1, j+self._params["max_shift"]+1) ):
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
    def decode_trace(self,trace=None,show_structures=False):
        if trace is None:
            trace = self.traceback()

        rnas = (self.rnaA,self.rnaB,self.rnaB)
        pos = [0]*3
        alignment = [""]*3
        for i,y in enumerate(trace):
            for s in range(3):
                if (y[s]==0):
                    alignment[s] = alignment[s] + "-"
                elif (y[s]==1):
                    alignment[s] = alignment[s] + rnas[s]["seq"][pos[s]]
                    pos[s]+=1

        if not show_structures:
            return alignment

        # annotate with structure
        anno_alignment = list()
        for alistr,rna in zip(alignment,rnas):
            anno_alignment.append( self._transfer_gaps(alistr,rna["structure"]) )
            anno_alignment.append( alistr )

        return anno_alignment

    # evaluate trace
    def eval_trace(self,trace=None):
        if trace is None:
            trace = self.traceback()
        pos=[0]*3
        for i,y in enumerate(trace):
            for k in range(3): pos[k] += y[k]
            # lookup case
            for x in self.recursionCases(pos[0],pos[1],pos[2]):
                if x[0] == y:
                    line = " ".join([ str(x) for x in 
                            [pos,
                            y,
                            x[1],
                            "-->",
                            self.evalCase(x,pos[0],pos[1],pos[2])
                            ]
                        ])
                    yield line
                    break

#convenience function to construct bialigner, perform and return text output
def bialign_triplet(seqA,seqB,strA,strB,show_structures,verbose,**args):
    ba = BiAlignerTriplet(seqA,seqB,strA,strB,**args)

    optscore = ba.optimize()
    yield "SCORE:" + str(optscore)

    for s in ba.decode_trace(show_structures=show_structures):
        yield s

    if verbose:
        yield from ba.eval_trace()

def main():
    parser = argparse.ArgumentParser(description= "Bialignment.")
    add_bialign_parameters(parser)
    
    args = parser.parse_args()
    print("Args",vars(args))

    for line in bialign_triplet(**vars(args)):
        print(line)

if __name__ == "__main__":
    main()
