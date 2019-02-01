#!/usr/bin/env python3
import RNA

import argparse
import numpy as np
from math import log,exp

## Alignment factory
class BiAligner:
   def __init__(self, seqA, seqB):
      self.rnaA = self._predict(seqA)
      self.rnaB = self._predict(seqB)

      # some parametrization
      self._sequence_match_similarity = 100
      self._sequence_mismatch_similarity = -50
      self._structure_weight = 200
      self._gap_cost = -50
      self._shift_cost = -50

      # precompute expected pairing partner for structure scores
      self.rnaA["ep"] = self._expected_pairing(self.rnaA)
      self.rnaB["ep"] = self._expected_pairing(self.rnaB)

      #print("RNA_A",self.rnaA)
      #print("RNA_B",self.rnaB)

      # the dynamic programming matrix
      self.M = None

   # iterator over recursion cases
   # 
   # per case:
   #       pair of access info
   #     and 
   #       function that returns list of case score components
   def recursionCases(self,i,j,k):
      # synchronous cases
      yield ((1,1,1), [self.mu1(i,j), self.mu2(i,k)])
      yield ((1,0,0), [self.g1A(i),   self.g2A(i)])
      yield ((0,1,1), [self.g1B(j),   self.g2B(k)])
      # shifting
      yield ((1,1,0), [self.mu1(i,j), self.g2A(i), self.Delta()])
      yield ((1,0,1), [self.mu2(i,k), self.g1A(i), self.Delta()])
      yield ((0,1,0), [self.g1A(i),   self.Delta()])
      yield ((0,0,1), [self.g2B(k),   self.Delta()])

   # plus operator (max in optimization; sum in pf)
   def plus(self, xs):
      if xs != []:
         return max(xs)
      else:
         return 0

   # mul operator (sum in optimization, product in pf)
   def mul(self, xs):
      return sum(xs)

   def guardCase(self,x,i,j,k):
      (io,jo,ko) = x[0]
      return i-io>=0 and j-jo>=0 and k-ko>=0

   def evalCase(self,x,i,j,k):
      (io,jo,ko) = x[0]
      return self.mul([self.M[i-io,j-jo,k-ko],
                       self.mul(x[1])
                      ])

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
   def _predict(sequence):
      prediction = dict()
      prediction["seq"] = str(sequence)
      prediction["len"] = len(prediction["seq"])
      
      fc = RNA.fold_compound(str(sequence))
      prediction["mfe"] = fc.mfe()
      prediction["pf"] = fc.pf()
      prediction["sbpp"] = BiAligner._symmetrize_bpps( fc.bpp() )
      return prediction

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
                  exp(-0.5+1/(1+exp(abs( self.rnaA["ep"][i] - self.rnaB["ep"][j] )))))
     
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

   # cost of shifting the 2 scores against each other
   def Delta(self):
      return self._shift_cost

   # run alignment algorithm
   def optimize(self):
      lenA = self.rnaA["len"]
      lenB = self.rnaB["len"]

      self.M = np.zeros((lenA+1,lenB+1,lenB+1), dtype=int)

      for i in range(0,lenA+1):
         for j in range(0,lenB+1):
            for k in range(0,lenB+1):
               self.M[i,j,k] = self.plus([ self.evalCase(x,i,j,k) 
                                           for x in self.recursionCases(i,j,k)
                                           if self.guardCase(x,i,j,k) ])
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
               print(y,
                     x[1],
                     "-->",
                     self.evalCase(x,pos[0],pos[1],pos[2]))
               break
         

def main(args):
   ba = BiAligner(args.seqA,args.seqB)

   optscore = ba.optimize()
   print("SCORE:",optscore)
   trace    = ba.traceback()
   print("TRACE:")
   for s in ba.decode_trace(trace): print("   "+s)
   ba.eval_trace(trace)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= "Bialignment.")
    parser.add_argument("seqA",help="RNA sequence A")
    parser.add_argument("seqB",help="RNA sequence B")
    
    args = parser.parse_args()
    main(args)
