#!/usr/bin/env python3

import RNA

import itertools
import sys
import argparse
import numpy as np
from math import log,exp,sqrt

VERSION_STRING = "BiAlign 0.3a"


class SparseMatrix4D:
    """Sparse 4D matrix of integers

    Valid entries (i,j,k,l) satisfy
    * i in range(n+1)
    * j in range(m+1)
    * k in range(i-max_shift, i+max_shift+1)
    * l in range(j-max_shift, j+max_shift+1)
    """

    def __init__(self, n, m, max_shift):
        self._n, self._m, self._max_shift = n, m, max_shift
        self._M = np.zeros((self._n+1, self._m+1,
            2*self._max_shift+1, 2*self._max_shift+1),
            dtype=int)

    def index(self, k):
        k = list(k)
        k[2] -= k[0] - self._max_shift
        k[3] -= k[1] - self._max_shift
        return tuple(k)

    def __getitem__(self, key):
        return self._M[self.index(key)]

    def __setitem__(self,key,value):
        self._M[self.index(key)] = value

class AffineDPMatrices:
    def __init__(self, n, m, max_shift):
        self._states = [ key for key in itertools.product(range(2),repeat=4)
            if (key[0]!=0 or key[1]!=0) and (key[2]!=0 or key[3]!=0)]
        self._Ms = {key:SparseMatrix4D(n, m, max_shift) for key in self._states}

    def __getitem__(self, key):
        return self._Ms[key]

    @property
    def states(self):
        return self._states

## Alignment factory
class BiAligner:
    nl = 14
    outmodes = {
        "sorted": [0, 1, 5, 3, 2, 4, nl] + [7, 6, 10, 8, 9, 11, nl] + [12, 13],
        "sorted_sym": [0, 1, 3, 2, 5, 4, nl] + [6, 7, 9, 8, 11, 10, nl] + [12, 13],
        "sorted_terse": [1, 5, 3, 4, nl] + [6, 10, 8, 11, nl] +
        [12, 13],
        "raw": [1,3,7,9],
        "raw_struct": list(range(4)) + list(range(6,10)),
        "full": range(nl)
        }

    def __init__(self, seqA, seqB, strA, strB, **params):
        # parametrization
        self._params = params

        self.molA = self._preprocess_seq(seqA, strA)
        self.molB = self._preprocess_seq(seqB, strB)

        self.gamma = self._params["gap_cost"]
        self.beta = self._params["gap_opening_cost"]

        if self._params["simmatrix"]:
            self._simmatrix = read_simmatrix(self._params["simmatrix"])
        else:
            self._simmatrix = None

        self._max_shift = self._params["max_shift"]


        #print( "Structure similarity",
        #       [ ((i,j),s) for i in range(1,len(seqA)+1)
        #                   for j in range(1,len(seqB)+1)
        #                   for s in [self._structure_similarity(i,j)] if s!=0])

        # the dynamic programming matrix
        self._M = None

    @property
    def _is_rna(self):
        return self._params['type'] == "RNA"

    @property
    def _affine(self):
        return self.beta != 0

    @staticmethod
    def error(text):
        print('ERROR:', text)
        sys.exit(-1)

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
    def recursion_cases(self,i,j,k,l):
        """recursion cases for non-affine alignment"""
        mu1ij = self.mu1(i,j)
        mu2kl = self.mu2(k,l)
        Delta = self._params["shift_cost"]

        # synchronous cases
        yield ((1,1,1,1), mu1ij + mu2kl)
        yield ((1,0,1,0), self.gamma   + self.gamma)
        yield ((0,1,0,1), self.gamma   + self.gamma)
        # shifting
        yield ((1,1,0,0), mu1ij + Delta)
        yield ((0,0,1,1), mu2kl + Delta)

        yield ((1,0,0,0), self.gamma + Delta)
        yield ((0,1,0,0), self.gamma + Delta)
        yield ((0,0,1,0), self.gamma + Delta)
        yield ((0,0,0,1), self.gamma + Delta)

        yield ((1,0,1,1), self.gamma + mu2kl + Delta)
        yield ((0,1,1,1), self.gamma + mu2kl + Delta)
        yield ((1,1,1,0), self.gamma + mu1ij + Delta)
        yield ((1,1,0,1), self.gamma + mu1ij + Delta)

        # double-shift cases -- these cases can be replaced by two others --> skip
        # yield ((0,1,1,0), self.gamma + self.gamma + 2 * self._params["shift_cost"])
        # yield ((1,0,0,1), self.gamma + self.gamma + 2 * self._params["shift_cost"])

    def affine_recursion_cases(self, state, idx):
        """yields recursion cases and their cost for affine gap cost
        """
        i, j, k, l = idx

        Delta = self._params["shift_cost"]
        mu1 = self.mu1(i,j)
        mu2 = self.mu2(k,l)

        def cost(source_state, x):
            cost = Delta * (abs(x[0] - x[2]) + abs(x[1] - x[3]))

            for a, b, mu in [(0, 1, mu1), (2, 3, mu2)]:
                if x[a] and x[b]: # match
                    cost += mu
                elif not x[a] and not x[b]:
                    pass
                elif x[a]:
                    cost += self.gamma
                    if source_state[a:b+1] != (1,0):
                        cost += self.beta # gap opening
                elif x[b]:
                    cost += self.gamma
                    if source_state[a:b+1] != (0,1):
                        cost += self.beta # gap opening
            return cost

        if self.guard_case(state, idx):
            for source_state in self.states:
                yield (source_state, state, cost(source_state, state))

        half_states = [(1,1), (1,0), (0,1)]

        offset = (0, 0, state[2], state[3])
        if self.guard_case(offset,idx):
            for half_state in half_states:
                target_state = (state[0], state[1], half_state[0], half_state[1])
                yield (target_state, offset, cost(target_state, offset))
        offset = (state[0], state[1], 0, 0)
        if self.guard_case(offset,idx):
            for half_state in half_states:
                target_state = (half_state[0], half_state[1], state[2], state[3])
                yield (target_state, offset, cost(target_state, offset))

    # plus operator (max in optimization; sum in pf)
    def plus(self, xs):
        try:
            return max(xs)
        except ValueError:
            return -1<<30

    ## mul operator (sum in optimization, product in pf)
    # def mul(self, xs):
    #    return sum(xs)

    def guard_case(self, offset, idx):
        (i,j,k,l) = idx
        (io,jo,ko,lo) = offset
        return i-io>=0 and j-jo>=0 and k-ko>=0 and l-lo>=0 and abs(k-ko-(i-io))<=self._max_shift and abs(l-lo-(j-jo))<=self._max_shift

    def eval_case(self, x, idx):
        i, j, k, l = idx
        io, jo, ko, lo = x[0]
        return self._M[i-io,j-jo,k-ko,l-lo] + x[1]

    def affine_eval_case(self, x, idx):
        i, j, k, l = idx
        io, jo, ko, lo = x[1]
        return self._M[x[0]][i-io,j-jo,k-ko,l-lo] + x[2]

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

    def _preprocess_seq(self, sequence, structure):
        x = dict()
        x["seq"] = str(sequence)
        x["len"] = len(x["seq"])

        if structure is None:
            if self._is_rna:
                fc = RNA.fold_compound(str(sequence))
                x["mfe"] = fc.mfe()
                x["pf"] = fc.pf()
                x["sbpp"] = BiAligner._symmetrize_bpps( fc.bpp() )
                x["mea"] = mea(x["sbpp"])
                x["structure"] = x["pf"][0]
            else:
                self.error("Structures have to be provided when aligning proteins")
        else:
            if len(structure)!=len(sequence):
                self.error("Provided structure and sequence must have the same length.")
            x["structure"] = structure
            if self._is_rna:
                x["sbpp"] = BiAligner._bp_matrix_from_fixed_structure(structure)

        n = x["len"]

        #note: the arrays are 1-based (we just ignore entries at 0)
        if self._is_rna:
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
    def _expected_pairing(mol):
        n = mol["len"]
        sbpp = mol["sbpp"]
        def ep(i):
            return sum( sbpp[i,j]*(j-i) for j in range(1,n+1) )

        return [0] + [ ep(i) for i in range(1,n+1) ]

    # sequence similarity of residues i and j, 1-based
    def _sequence_similarity(self,i,j):
        if self._simmatrix:
            return self._simmatrix[self.molA["seq"][i-1]][self.molB["seq"][j-1]]

        if self.molA["seq"][i-1]==self.molB["seq"][j-1]:
                return self._params["sequence_match_similarity"]
        else:
            return self._params["sequence_mismatch_similarity"]

    def _structure_similarity(self,i,j):
        if self._is_rna:
            sim = int( self._params["structure_weight"] *
                        (
                          sqrt(self.molA["up"][i]*self.molB["up"][j])
                          + sqrt(self.molA["down"][i]*self.molB["down"][j])
                          + sqrt(self.molA["unp"][i]*self.molB["unp"][j])
                        )
                      )
        else:
            if self.molA["structure"][i-1]==self.molB["structure"][j-1]:
                sim = 100
            else:
                sim = 0
        return sim

    # Scoring functions
    # note: scoring functions have 1-based indices

    # match/mismatch cost for i~j, score 1
    def mu1(self, i, j):
        return self._sequence_similarity(i, j)

    # match/mismatch cost for i~j, score 2
    def mu2(self, i, j):
        return self._structure_similarity(i, j)

    # run non-affine alignment algorithm
    def optimize(self):
        if self._affine:
            return self.affine_optimize()

        lenA = self.molA["len"]
        lenB = self.molB["len"]

        self._M = np.zeros((lenA+1,lenB+1,lenA+1,lenB+1), dtype='int32')

        for i in range(0,lenA+1):
            for j in range(0,lenB+1):
                for k in range( max(0, i-self._max_shift), min(lenA+1, i+self._max_shift+1) ):
                    for l in range( max(0, j-self._max_shift), min(lenB+1, j+self._max_shift+1) ):
                        idx = (i, j, k, l)
                        if idx == (0,0,0,0):
                            continue
                        self._M[idx] = self.plus(
                            self.eval_case(x, idx)
                            for x in self.recursion_cases(idx)
                            if self.guard_case(x[0],idx))
        return self._M[lenA,lenB,lenA,lenB]

    # run affine alignment algorithm
    def affine_optimize(self):
        lenA = self.molA["len"]
        lenB = self.molB["len"]

        self._M = AffineDPMatrices(lenA, lenB, self._max_shift)
        self.states = self._M.states

        # initialize - [0,0,0,0] is finite only if state is 'both match'
        for state in self.states:
            self._M[state][0,0,0,0] = -1<<30 # -infinity
        self._M[(1,1,1,1)][0,0,0,0] = 0

        for i in range(0,lenA+1):
            for j in range(0,lenB+1):
                for k in range( max(0, i-self._max_shift), min(lenA+1, i+self._max_shift+1) ):
                    for l in range( max(0, j-self._max_shift), min(lenB+1, j+self._max_shift+1) ):
                        idx = (i,j,k,l)
                        if idx == (0,0,0,0): # "initialization"
                            continue
                        for state in self.states:
                           self._M[state][idx] = self.plus(
                               self.affine_eval_case(x, idx)
                               for x in self.affine_recursion_cases(state, idx))

        return max(self._M[state][lenA,lenB,lenA,lenB] for state in self.states)

    # perform traceback
    # @returns list of 'trace arrows'
    def traceback(self):
        if self._affine:
            return self.affine_traceback()
        lenA = self.molA["len"]
        lenB = self.molB["len"]

        trace=[]

        def trace_from(i, j , k, l):
            for x in self.recursion_cases(i,j,k,l):
                if self.guard_case(x[0], (i,j,k,l)):
                    if self.eval_case(x, (i, j, k, l)) == self._M[i,j,k,l]:
                        (io,jo,ko,lo) = x[0]
                        trace.append((io,jo,ko,lo))
                        trace_from(i-io,j-jo,k-ko,l-lo)
                        break

        trace_from(lenA,lenB,lenA,lenB)
        return list(reversed(trace))

    # perform traceback
    # @returns list of 'trace arrows'
    def affine_traceback(self):
        lenA = self.molA["len"]
        lenB = self.molB["len"]

        trace=[]

        def trace_from(state, idx):
            if idx == (0,0,0,0) and state == (1,1,1,1):
                return True
            i, j, k, l = idx
            for x in self.affine_recursion_cases(state, idx):
                if self.guard_case(x[1], idx):
                    if self.affine_eval_case(x, idx) == self._M[state][idx]:
                        (io,jo,ko,lo) = x[1]
                        trace.append((io, jo, ko, lo))
                        return trace_from(x[0], (i-io,j-jo,k-ko,l-lo))
            return False
        myidx = np.argmax([self._M[state][lenA,lenB,lenA,lenB] for state in self.states])
        best_state=self.states[myidx]

        if not trace_from(best_state, (lenA, lenB, lenA, lenB)):
            print("WARNING: incomplete traceback. Alignment could be garbage.")
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

        mols = (self.molA,self.molB,self.molA,self.molB)
        pos = [0]*len(mols)
        alignment = [""]*len(mols)
        for i,y in enumerate(trace):
            for s in range(len(mols)):
                if (y[s]==0):
                    alignment[s] = alignment[s] + "-"
                elif (y[s]==1):
                    alignment[s] = alignment[s] + mols[s]["seq"][pos[s]]
                    pos[s]+=1

        #alignment[0],alignment[1] = highlight_sequence_identity(alignment[0],alignment[1])
        #alignment[2],alignment[3] = highlight_sequence_identity(alignment[2],alignment[3])

        # compute consensus sequences
        cons_seq = [consensus_sequence(alignment[2*i],alignment[2*i+1])
                    for i in range(2)]

        # annotate with structure
        anno_ali = list()
        for alistr,mol in zip(alignment,mols):
            anno_ali.append( self._transfer_gaps(alistr,mol["structure"]) )
            anno_ali.append( alistr )

        for i,j in [(4,6),(0,2)]:
            if self._is_rna:
                sbpp = consensus_sbpp( alistrA = anno_ali[i], alistrB = anno_ali[j],
                                       sbppA=self.molA["sbpp"], sbppB=self.molB["sbpp"]                                     )
                structure = mea(sbpp, brackets="[]")[0]
            else:
                structure = consensus_sequence(anno_ali[i], anno_ali[j])
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

        mode = self.auto_complete(self._params["outmode"],self.outmodes.keys())

        if mode in self.outmodes:
            order = self.outmodes[mode]
        else:
            print("WARNING: unknown output mode. Expect one of " +
                    str(list(self.outmodes.keys())) )
            order = self.outmodes["sorted"]

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
            for x in self.recursion_cases(pos[0],pos[1],pos[2],pos[3]):
                if x[0] == y:
                    line = " ".join([ str(x) for x in
                            [pos,
                            y,
                            x[1],
                            "-->",
                            self.eval_case(x, (pos[0],pos[1],pos[2],pos[3]))
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
    parser.add_argument("seqA",help="sequence A")
    parser.add_argument("seqB",help="sequence B")
    parser.add_argument("--strA",default=None,help="structure A")
    parser.add_argument("--strB",default=None,help="structure B")
    parser.add_argument("--nameA",default="A",help="name A")
    parser.add_argument("--nameB",default="B",help="name B")
    parser.add_argument("-v","--verbose",action='store_true',help="Verbose")

    parser.add_argument("--type", default="RNA", type=str,
        help="Type of molecule: RNA or Protein")

    parser.add_argument("--nodescription",action='store_true',
                        help="Don't prefix the strings in output alignment with descriptions")
    parser.add_argument("--outmode", default="sorted",
                        help="Output mode [call --mode help for a list of options]")

    parser.add_argument("--sequence_match_similarity", type=int,
            default=100, help="Similarity of matching nucleotides")
    parser.add_argument("--sequence_mismatch_similarity", type=int,
            default=0, help="Similarity of mismatching nucleotides")
    parser.add_argument("--structure_weight", type=int, default=100,
            help="Weighting factor for structure similarity")
    parser.add_argument("--gap_opening_cost", type=int, default=0,
            help="Similarity of opening a gap (turns on affine gap cost if not 0)")
    parser.add_argument("--gap_cost", type=int, default=-200,
            help="Similarity of a single gap position")
    parser.add_argument("--shift_cost", type=int, default=-250,
            help="Similarity of shifting the two scores against each other")
    parser.add_argument("--max_shift", type=int, default=2,
            help="Maximal number of shifts away from the diagonal in either direction")

    parser.add_argument("--version", action='version', version=VERSION_STRING )

    parser.add_argument("--simmatrix", type=str, default = None,
            help = "Similarity matrix")

blosum62 =\
"""-  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

def read_simmatrix(filename, scale = 100):
    if filename == "BLOSUM62":
        lines = blosum62.split('\n')
    else:
        with open(filename,'r') as fh:
            lines = fh.readlines()

    keys = None
    keys2 = []
    matrix = dict()

    for i, line in enumerate(lines):
        if keys and i>len(keys):
            break
        line = line.split()
        if line[0]=='-':
            keys = line[1:]
        else:
            keys2.append(line[0])
            matrix[line[0]] = {key:(scale * int(val)) for key, val in zip(keys,line[1:])}

    if not keys == keys2:
        print("ERROR while reading simmatrix {filename}.")
    return matrix

def main():
    parser = argparse.ArgumentParser(description= "Bialignment.")
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

    print("\n".join(["Input:","seqA\t "+args.seqA,"seqB\t "+args.seqB,"strA\t "+args.strA,"strB\t "+args.strB,""]))

    if args.outmode == "help":
        print()
        print("Available modes: "+", ".join(BiAligner.outmodes.keys()))
        print()
        exit()

    for line in bialign(**vars(args)):
        print(line)

if __name__ == "__main__":
    main()
