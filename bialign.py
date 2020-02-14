#!/usr/bin/env python3

import RNA

import sys
import argparse
import numpy as np
from math import log,exp,sqrt


## Alignment factory
class BiAligner:
    nl = 14
    modes = {
             "sorted": [0, 1, 5, 3, 2, 4, nl] + [7, 6, 10, 8, 9, 11, nl] + [12, 13],
             "sorted_sym": [0, 1, 3, 2, 5, 4, nl] + [6, 7, 9, 8, 11, 10, nl] + [12, 13],
             "sorted_terse": [1, 5, 3, 4, nl] + [6, 10, 8, 11, nl] +
             [12, 13],
             "raw": [1,3,7,9],
             "raw_struct": list(range(4)) + list(range(6,10)),
             "full": range(nl)
            }


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
        mu1ij = self.mu1(i,j)
        mu2kl = self.mu2(k,l)
        Delta = self._params["shift_cost"]

        # synchronous cases
        yield ((1,1,1,1), mu1ij + mu2kl)
        yield ((1,0,1,0), self.g1A(i)   + self.g2A(k))
        yield ((0,1,0,1), self.g1B(j)   + self.g2B(l))
        # shifting
        yield ((1,1,0,0), mu1ij + Delta)
        yield ((0,0,1,1), mu2kl + Delta)

        yield ((1,0,0,0), self.g1A(i) + Delta)
        yield ((0,1,0,0), self.g1B(j) + Delta)
        yield ((0,0,1,0), self.g2A(k) + Delta)
        yield ((0,0,0,1), self.g2B(l) + Delta)

        yield ((1,0,1,1), self.g1A(i) + mu2kl + Delta)
        yield ((0,1,1,1), self.g1B(j) + mu2kl + Delta)
        yield ((1,1,1,0), self.g2A(k) + mu1ij + Delta)
        yield ((1,1,0,1), self.g2B(l) + mu1ij + Delta)

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
    #
    # NOTE: bpp and sbpp have 1-based access (row and column 0 are ignored) 
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
            x["mea"] = mea(x["sbpp"])
            x["structure"] = x["pf"][0]
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
        sim = int( self._params["structure_weight"] *
                    (
                      sqrt(self.rnaA["up"][i]*self.rnaB["up"][j])
                      + sqrt(self.rnaA["down"][i]*self.rnaB["down"][j])
                      + sqrt(self.rnaA["unp"][i]*self.rnaB["unp"][j])
                    )
                  )
        return sim
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

    @staticmethod
    def _shift_string(ali, idx):
        length = len(ali[0])
        
        def shift(i):
            c1 = "X"
            c2 = "X"
            if ali[idx][i]=="-":
                c1 = "-"
            if ali[idx+2][i]=="-":
                c2 = "-"

            if c1==c2:
                return "."
            elif c1=="-":
                return ">"
            elif c2=="-":
                return "<"

        s=[shift(i) for i in range(length)]
        return "".join(s)  

    @staticmethod
    def auto_complete(x,xs):
        xs = list(xs)
        xs.sort()
        for y in xs:
            if y.startswith(x):
                return y
        return x

    # decode trace to alignment strings
    def decode_trace(self, *, trace=None):
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

        #alignment[0],alignment[1] = highlight_sequence_identity(alignment[0],alignment[1])
        #alignment[2],alignment[3] = highlight_sequence_identity(alignment[2],alignment[3])

        # compute consensus sequences
        cons_seq = [consensus_sequence(alignment[2*i],alignment[2*i+1])
                    for i in range(2)]

        # annotate with structure
        anno_ali = list()
        for alistr,rna in zip(alignment,rnas):
            anno_ali.append( self._transfer_gaps(alistr,rna["structure"]) )
            anno_ali.append( alistr )

        for i,j in [(4,6),(0,2)]:
            sbpp = consensus_sbpp( alistrA = anno_ali[i], alistrB = anno_ali[j],
                                   sbppA=self.rnaA["sbpp"], sbppB=self.rnaB["sbpp"]                                     )
            structure = mea(sbpp, brackets="[]")[0]
            anno_ali.insert(j+2,structure)

        shift_strings=list()
        for i in range(2):
            shift_strings.append(self._shift_string(alignment,i))

        alignment = anno_ali

        alignment.insert(len(alignment),cons_seq[1])
        alignment.insert(len(alignment)//2,cons_seq[0])

        alignment.extend(shift_strings)

        nameA = self._params["nameA"]
        nameB = self._params["nameB"]

        struct_postfix=" ss"
        names = [
         nameA+struct_postfix,
         nameA,
         nameB+struct_postfix,
         nameB,
         "consensus"+struct_postfix,
         "consensus",
         nameA+struct_postfix,
         nameA,
         nameB+struct_postfix,
         nameB,
         "consensus"+struct_postfix,
         "consensus",
         nameA+" shifts",
         nameB+" shifts"
        ]

        width = max(map(len,names)) + 4

        #add names of single lines
        if not self._params["nodescription"]:
            alignment = [("{:{width}}{}").format(name, alistr, width=width) for alistr, name in zip(alignment, names)]

        alignment.append("")

        mode = self.auto_complete(self._params["mode"],self.modes.keys())

        if mode in self.modes:
            order = self.modes[mode]
        else:
            print("WARNING: unknown output mode. Expect one of " +
                    str(list(self.modes.keys())) )
            order = self.modes["sorted"]

        # re-sort
        alignment = [alignment[i] for i in order]
                    
        # re-sort, terse
        #alignment = [alignment[i] for i in ]


        return alignment

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

# compute mea structure
def mea(sbpp,gamma=3,*,brackets="()"):
    n = len(sbpp)-1
    
    F = np.zeros((n+1,n+1),dtype='float')
    T = np.zeros((n+1,n+1),dtype='int')

    candidates = [ list() for i in range(0,n+1) ]

    for i in reversed(range(1,n+1)):
        candidates[i].append((i, sbpp[i,i]))
        F[i,i-1]=0
        for j in range(i,n+1):
            for k,C in candidates[j]:
                if F[i,j] < F[i,k-1] + C:
                    F[i,j] = F[i,k-1] + C
                    T[i,j] = k

            if i+3>=j:
                continue

            C = F[i+1,j-1] + 2*gamma*sbpp[i,j]
            if C > F[i,j]:
                candidates[j].append((i,C))
                F[i,j] = C
                T[i,j] = i

    # trace back
    
    structure = ['.']*(n+1)
    stack = list()
    stack.append((1,n))
    
    while stack:
        (i,j) = stack.pop()
        k = T[i,j]
        if i+3>=j or k==0: continue

        if k == j:
            stack.append((i,j-1))
        elif k == i:
            structure[k]=brackets[0]
            structure[j]=brackets[1]
            stack.append((k+1,j-1))
        else:
            stack.append((i,k-1))
            stack.append((k+1,j-1))
            structure[k]=brackets[0]
            structure[j]=brackets[1]

    return ("".join(structure[1:]),F[1,n])

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

def consensus_sequence(alistrA,alistrB):
    def c(x,y):
        if (x==y):
            return x # "*"
        else:
            return "."
    return "".join([c(x,y) for x,y in
                   zip(alistrA.upper(),alistrB.upper())])

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

# consensus base pair probabilities
def consensus_sbpp(alistrA,sbppA,alistrB,sbppB):
    sbpp = np.zeros( (len(alistrA)+1, len(alistrB)+1), dtype='float' )
    
    length = [len(sbppA)-1, len(sbppB)-1]

    p0 = [1,1]
    for i0,x0 in enumerate(zip(alistrA, alistrB)):
        p1 = [1,1]
        for i1,x1 in enumerate(zip(alistrA, alistrB)):
            pr = [0,0]
            for k,sbppX in [(0,sbppA),(1,sbppB)]:
                if x0[k]=='-' or x1[k]=='-':
                    pr[k]=0
                else:
                    pr[k] = sbppX[p0[k],p1[k]]

            sbpp[i0+1,i1+1] = sqrt( pr[0] * pr[1] )
            for k in range(2):
                if x1[k]!='-': 
                    p1[k]+=1
        for k in range(2):
            if x0[k]!='-':
                p0[k]+=1

    return sbpp

# highlight matched base pairs in a pairwise alignment; for balanced dot bracket strings
def highlight_structure_identity(alistrA, alistrB):

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

# highlight matched base pairs in a pairwise alignment; given sbpp matrices
def highlight_structure_similarity(alistrA, alistrB, *, sbppA, sbppB):

    sbpp = consensus_sbpp(alistrA,sbppA,alistrB,sbppB)

    structure = parse_dotbracket(mea(sbpp)[0])

    res = [ list(x) for x in [alistrA, alistrB] ]
    for i in range(len(alistrA)):
        for j in range(i+1,len(alistrA)): 
            if structure[i] == j: 
                res[0][i]="<"
                res[1][i]="<"
                res[0][j]=">"
                res[1][j]=">"

    return [ "".join(x) for x in res]

def bialign(seqA, seqB, strA, strB, verbose, **args):
    ba = BiAligner(seqA,seqB,strA,strB,**args)

    optscore = ba.optimize()
    yield "SCORE: " + str(optscore)
    yield ""

    ali = ba.decode_trace()

    yield from ali

    if verbose:
        yield from ba.eval_trace()

def add_bialign_parameters(parser):
    parser.add_argument("seqA",help="RNA sequence A")
    parser.add_argument("seqB",help="RNA sequence B")
    parser.add_argument("--strA",default=None,help="RNA structure A")
    parser.add_argument("--strB",default=None,help="RNA structure B")
    parser.add_argument("--nameA",default="RNA A",help="RNA name A")
    parser.add_argument("--nameB",default="RNA B",help="RNA name B")
    parser.add_argument("-v","--verbose",action='store_true',help="Verbose")

    parser.add_argument("--nodescription",action='store_true',
                        help="Don't prefix the strings in output alignment with descriptions")
    parser.add_argument("--mode", default="sorted",
                        help="Output mode [call --mode help for a list of options]")

    parser.add_argument("--sequence_match_similarity", type=int,
            default=100, help="Similarity of matching nucleotides")
    parser.add_argument("--sequence_mismatch_similarity", type=int,
            default=0, help="Similarity of mismatching nucleotides")
    parser.add_argument("--structure_weight", type=int, default=100,
            help="Weighting factor for structure similarity")
    parser.add_argument("--gap_cost", type=int, default=-200,
            help="Similarity of a single gap")
    parser.add_argument("--shift_cost", type=int, default=-250,
            help="Similarity of shifting the two scores against each other")
    parser.add_argument("--max_shift", type=int, default=2,
            help="Maximal number of shifts away from the diagonal in either direction")

def main():
    parser = argparse.ArgumentParser(description= "Bialignment.")
    add_bialign_parameters(parser)

    args = parser.parse_args()

    if args.mode == "help":
        print()
        print("Available modes: "+", ".join(BiAligner.modes.keys()))
        print()
        exit()

    for line in bialign(**vars(args)):
        print(line)

if __name__ == "__main__":
    main()
