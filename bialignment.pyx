#!/usr/bin/env python3


from bialignment_nonpyx import *
import itertools
import sys
import numpy as np
from math import log, exp, sqrt

__version__ = "0.3a"


cdef class SparseMatrix4D:
    """Sparse 4D matrix of integers

    Valid entries (i,j,k,l) satisfy
    * i in range(n+1)
    * j in range(m+1)
    * k in range(i-max_shift, i+max_shift+1)
    * l in range(j-max_shift, j+max_shift+1)
    """

    cdef int _n, _m, max_shift
    cdef _M

    def __init__(self, n, m, max_shift):
        self._n, self._m, self.max_shift = n, m, max_shift
        self._M = np.zeros(
            (
                self._n + 1,
                self._m + 1,
                2 * self.max_shift + 1,
                2 * self.max_shift + 1,
            ),
            dtype=int,
        )

    def __getitem__(self, k):
        return self._M[k[0], k[1], k[2]-k[0]+self.max_shift, k[3]-k[1]+self.max_shift]

    def __setitem__(self, k, value):
        self._M[k[0], k[1], k[2]-k[0]+self.max_shift, k[3]-k[1]+self.max_shift] = value


cdef fourtuple_minus(int x[4], int y[4], int z[4]):
    for i in range(4):
      z[i] = x[i] - y[i]

cdef int fourbits_to_int(int key[4]):
    cdef int i
    cdef int res = 0
    for i in range(4):
        res = res << 1
        res += key[i]
    return res

cdef class AffineDPMatrices:
    cdef _Ms
    cdef _states

    def __init__(self, n, m, max_shift):
        self._states = [
            key
            for key in itertools.product(range(2), repeat=4)
            if (key[0] != 0 or key[1] != 0) and (key[2] != 0 or key[3] != 0)
        ]
        self._Ms = [None]*16
        cdef int ckey[4]
        for key in self._states:
            ckey = key
            idx = fourbits_to_int(ckey)
            self._Ms[idx] = SparseMatrix4D(n, m, max_shift)

    def __getitem__(self, key):
        cdef int ckey[4]
        ckey = key
        idx = fourbits_to_int(ckey)
        return self._Ms[idx]

    @property
    def states(self):
        return self._states

## BiAligner affine score
cdef int affine_score(int source_state[4], int x[4], int mu1, int mu2, int beta, int gamma, int Delta):
    """
    Score of one alignment column in the affine case
    @param source_state state of prefix alignment without new column
    @param x new column as 0,1 vector (0=gap, 1=non-gap)
    @param mu1 similarity in case of seq match    
    @param mu2 similarity in case of str match
    @param beta gap opening
    @param gamma gap extension
    @param Delta shift cost
    
    @return score
    """
    cdef int score = Delta * (abs(x[0] - x[2]) + abs(x[1] - x[3]))

    cdef int a
    cdef int b
    cdef int mu

    a, b, mu = (0, 1, mu1)
    xa = x[a]
    xb = x[b]
    if xa and xb:  # match
        score += mu
    elif xa and not xb:
        score += gamma
        if not (source_state[a] == 1 and source_state[b] == 0):
            score += beta  # gap opening
    elif not xa and xb:
        score += gamma
        if not (source_state[a] == 0 and source_state[b] == 1):
            score += beta  # gap opening

    a, b, mu = (2, 3, mu2)
    xa = x[a]
    xb = x[b]
    if xa and xb:  # match
        score += mu
    elif xa and not xb:
        score += gamma
        if not (source_state[a] == 1 and source_state[b] == 0):
            score += beta  # gap opening
    elif not xa and xb:
        score += gamma
        if not (source_state[a] == 0 and source_state[b] == 1):
            score += beta  # gap opening

    return score

cdef int cguard_case(int o[4], int x[4], int max_shift):
    return (
        x[0] - o[0] >= 0
        and x[1] - o[1] >= 0
        and x[2] - o[2] >= 0
        and x[3] - o[3] >= 0
        and abs(x[2] - o[2] - (x[0] - o[0])) <= max_shift
        and abs(x[3] - o[3] - (x[1] - o[1])) <= max_shift
    )

def guard_case(o, x, int max_shift):
    cdef int co[4]
    co = o
    cdef int cx[4]
    cx = x
    return cguard_case(co, cx, max_shift)


def argmin(xs):
    return min(enumerate(xs), key=lambda x:x[1])[0]

## Alignment factory
cdef class BiAligner:

    cdef int beta
    cdef int gamma
    cdef int max_shift
    cdef _params
    cdef molA
    cdef molB
    cdef _M
    cdef _simmatrix
    cdef states
    cdef int cstates[9][4]

    nl = 14
    outmodes = {
        "sorted": [0, 1, 5, 3, 2, 4, nl] + [7, 6, 10, 8, 9, 11, nl] + [12, 13],
        "sorted_sym": [0, 1, 3, 2, 5, 4, nl] + [6, 7, 9, 8, 11, 10, nl] + [12, 13],
        "sorted_terse": [1, 5, 3, 4, nl] + [6, 10, 8, 11, nl] + [12, 13],
        "raw": [1, 3, 7, 9],
        "raw_struct": list(range(4)) + list(range(6, 10)),
        "full": range(nl),
    }

    def __init__(self, seqA, seqB, strA, strB, **params):
        # parametrization
        self._params = params

        self.molA = self._preprocess_seq(seqA, strA)
        self.molB = self._preprocess_seq(seqB, strB)

        self.gamma = self._params["gap_cost"]
        self.beta = self._params["gap_opening_cost"]
        self.max_shift = self._params["max_shift"]

        if self._params["simmatrix"]:
            self._simmatrix = read_simmatrix(self._params["simmatrix"])
        else:
            self._simmatrix = None


        # the dynamic programming matrix
        self._M = None

    @property
    def _is_rna(self):
        return self._params["type"] == "RNA"

    @property
    def _affine(self):
        return self.beta != 0

    @staticmethod
    def error(text):
        print("ERROR:", text)
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
    def recursion_cases(self, idx):
        """recursion cases for non-affine alignment"""
        i, j, k, l = idx
        mu1ij = self.mu1(i, j)
        mu2kl = self.mu2(k, l)
        Delta = self._params["shift_cost"]

        # synchronous cases
        yield ((1, 1, 1, 1), mu1ij + mu2kl)
        yield ((1, 0, 1, 0), self.gamma + self.gamma)
        yield ((0, 1, 0, 1), self.gamma + self.gamma)
        # shifting
        yield ((1, 1, 0, 0), mu1ij + Delta)
        yield ((0, 0, 1, 1), mu2kl + Delta)

        yield ((1, 0, 0, 0), self.gamma + Delta)
        yield ((0, 1, 0, 0), self.gamma + Delta)
        yield ((0, 0, 1, 0), self.gamma + Delta)
        yield ((0, 0, 0, 1), self.gamma + Delta)

        yield ((1, 0, 1, 1), self.gamma + mu2kl + Delta)
        yield ((0, 1, 1, 1), self.gamma + mu2kl + Delta)
        yield ((1, 1, 1, 0), self.gamma + mu1ij + Delta)
        yield ((1, 1, 0, 1), self.gamma + mu1ij + Delta)

        # double-shift cases -- these cases can be replaced by two others --> skip
        # yield ((0,1,1,0), self.gamma + self.gamma + 2 * self._params["shift_cost"])
        # yield ((1,0,0,1), self.gamma + self.gamma + 2 * self._params["shift_cost"])


    def affine_recursion_cases(self, state, idx):
        """yields recursion cases and their score for affine gap cost"""
        i, j, k, l = idx

        cdef int Delta = self._params["shift_cost"]
        cdef int mu1 = self.mu1(i, j)
        cdef int mu2 = self.mu2(k, l)
        cdef int beta = self.beta
        cdef int gamma = self.gamma

        cdef int csource_state[4]
        cdef int cstate[4]
        cstate = state
        cdef int cidx[4]
        cidx = idx

        cdef int max_shift = self.max_shift

        cdef int ss

        if cguard_case(cstate, cidx, max_shift):
            for ss in range(9):
                csource_state = self.cstates[ss]
                yield (csource_state, cstate,
                    affine_score(csource_state, cstate, mu1, mu2, beta, gamma, Delta))

        cdef int half_states[3][2]
        half_states = [[1, 1], [1, 0], [0, 1]]

        cdef int offset[4]
        offset = [0, 0, cstate[2], cstate[3]]
        if cguard_case(offset, cidx, max_shift):
            for hs in range(3):
                csource_state = (cstate[0], cstate[1], half_states[hs][0], half_states[hs][1])
                yield (csource_state, offset,
                    affine_score(csource_state, offset, mu1, mu2, beta, gamma, Delta))
        offset = [cstate[0], cstate[1], 0, 0]
        if cguard_case(offset, cidx, max_shift):
            for hs in range(3):
                csource_state = (half_states[hs][0], half_states[hs][1], cstate[2], cstate[3])
                yield (csource_state, offset,
                    affine_score(csource_state, offset, mu1, mu2, beta, gamma, Delta))

    # plus operator (max in optimization; sum in pf)
    def plus(self, xs):
        try:
            return max(xs)
        except ValueError:
            return -1 << 30

    ## mul operator (sum in optimization, product in pf)
    # def mul(self, xs):
    #    return sum(xs)


    def eval_case(self, x, idx):
        i, j, k, l = idx
        io, jo, ko, lo = x[0]
        return self._M[i - io, j - jo, k - ko, l - lo] + x[1]

    cdef int affine_eval_case(self, x, int idx[4]):
        cdef int x1[4]
        x1 = x[1]
        cdef int y[4]
        fourtuple_minus(idx, x1, y)
        return self._M[x[0]][y] + x[2]

    # make bpp symmetric (based on upper triangular matrix)
    # and set diagonal to unpaired probs
    #
    # NOTE: bpp and sbpp have 1-based access (row and column 0 are ignored)
    @staticmethod
    def _symmetrize_bpps(bpp):
        n = len(bpp) - 1
        sbpp = np.zeros((n + 1, n + 1), dtype="float")
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                sbpp[i, j] = bpp[i][j]
                sbpp[j, i] = bpp[i][j]

        for i in range(1, n + 1):
            sbpp[i, i] = 1.0 - sum(sbpp[i, j] for j in range(1, n + 1))

        return sbpp

    def _preprocess_seq(self, sequence, structure):
        x = dict()
        x["seq"] = str(sequence)
        x["len"] = len(x["seq"])

        if structure is None:
            if self._is_rna:
                import RNA
                fc = RNA.fold_compound(str(sequence))
                x["mfe"] = fc.mfe()
                x["pf"] = fc.pf()
                x["sbpp"] = BiAligner._symmetrize_bpps(fc.bpp())
                x["mea"] = mea(x["sbpp"])
                x["structure"] = x["pf"][0]
            else:
                self.error("Structures have to be provided when aligning proteins")
        else:
            if len(structure) != len(sequence):
                self.error("Provided structure and sequence must have the same length.")
            x["structure"] = structure
            if self._is_rna:
                x["sbpp"] = BiAligner._bp_matrix_from_fixed_structure(structure)

        n = x["len"]

        # note: the arrays are 1-based (we just ignore entries at 0)
        if self._is_rna:
            x["up"] = [
                sum(x["sbpp"][i][j] for j in range(1, i - 1)) for i in range(0, n + 1)
            ]
            x["down"] = [
                sum(x["sbpp"][i][j] for j in range(i + 1, n + 1))
                for i in range(0, n + 1)
            ]
            x["unp"] = [1.0 - x["up"][i] - x["down"][i] for i in range(0, n + 1)]

        return x

    @staticmethod
    def _bp_matrix_from_fixed_structure(structure):
        n = len(structure)
        bpm = np.zeros((n + 1, n + 1), dtype="float")
        stack = list()
        for i in range(n):
            if structure[i] == "(":
                stack.append(i)
            elif structure[i] == ")":
                j = stack.pop()
                bpm[i + 1, j + 1] = 1.0
                bpm[j + 1, i + 1] = 1.0
            else:
                bpm[i + 1, i + 1] = 1.0
        return bpm

    @staticmethod
    def _expected_pairing(mol):
        n = mol["len"]
        sbpp = mol["sbpp"]

        def ep(i):
            return sum(sbpp[i, j] * (j - i) for j in range(1, n + 1))

        return [0] + [ep(i) for i in range(1, n + 1)]

    # sequence similarity of residues i and j, 1-based
    def _sequence_similarity(self, i, j):
        if self._simmatrix:
            return self._simmatrix[self.molA["seq"][i - 1]][self.molB["seq"][j - 1]]

        if self.molA["seq"][i - 1] == self.molB["seq"][j - 1]:
            return self._params["sequence_match_similarity"]
        else:
            return self._params["sequence_mismatch_similarity"]

    def _structure_similarity(self, i, j):
        if self._is_rna:
            sim = int(
                self._params["structure_weight"]
                * (
                    sqrt(self.molA["up"][i] * self.molB["up"][j])
                    + sqrt(self.molA["down"][i] * self.molB["down"][j])
                    + sqrt(self.molA["unp"][i] * self.molB["unp"][j])
                )
            )
        else:
            if self.molA["structure"][i - 1] == self.molB["structure"][j - 1]:
                sim = self._params["structure_weight"]
            else:
                sim = 0
        return sim

    # Scoring functions
    # note: scoring functions have 1-based indices

    # match/mismatch score for i~j, score 1
    def mu1(self, i, j):
        return self._sequence_similarity(i, j)

    # match/mismatch score for i~j, score 2
    def mu2(self, i, j):
        return self._structure_similarity(i, j)

    # run non-affine alignment algorithm
    def optimize(self):
        if self._affine:
            return self.affine_optimize()

        cdef int max_shift = self.max_shift

        lenA = self.molA["len"]
        lenB = self.molB["len"]

        self._M = SparseMatrix4D(lenA, lenB, max_shift)

        for i in range(0, lenA + 1):
            for j in range(0, lenB + 1):
                for k in range(
                    max(0, i - self.max_shift), min(lenA + 1, i + self.max_shift + 1)
                ):
                    for l in range(
                        max(0, j - self.max_shift),
                        min(lenB + 1, j + self.max_shift + 1),
                    ):
                        idx = (i, j, k, l)
                        if idx == (0, 0, 0, 0):
                            continue
                        self._M[idx] = self.plus(
                            self.eval_case(x, idx)
                            for x in self.recursion_cases(idx)
                            if guard_case(x[0], idx, max_shift)
                        )
        return self._M[lenA, lenB, lenA, lenB]

    # run affine alignment algorithm
    def affine_optimize(self):
        cdef int lenA = self.molA["len"]
        cdef int lenB = self.molB["len"]

        self._M = AffineDPMatrices(lenA, lenB, self.max_shift)
        self.states = self._M.states
        self.cstates = self.states

        # initialize - [0,0,0,0] is finite only if state is 'both match'
        for state in self.states:
            self._M[state][0, 0, 0, 0] = -1 << 30  # -infinity
        self._M[(1, 1, 1, 1)][0, 0, 0, 0] = 0

        cdef int idx[4]
        cdef int i, j, k, l
        cdef int s

        for i in range(0, lenA + 1):
            for j in range(0, lenB + 1):
                for k in range(
                    max(0, i - self.max_shift), min(lenA + 1, i + self.max_shift + 1)
                ):
                    for l in range(
                        max(0, j - self.max_shift),
                        min(lenB + 1, j + self.max_shift + 1),
                    ):
                        idx = [i, j, k, l]
                        if idx == [0, 0, 0, 0]:  # "initialization"
                            continue
                        for s in range(9):# cstate in self.cstates:
                            self._M[self.cstates[s]][idx] = self.plus(
                                self.affine_eval_case(x, idx)
                                for x in self.affine_recursion_cases(self.cstates[s], idx)
                            )

        return max(self._M[state][lenA, lenB, lenA, lenB] for state in self.states)

    # perform traceback
    # @returns list of 'trace arrows'
    def traceback(self):
        if self._affine:
            return self.affine_traceback()
        lenA = self.molA["len"]
        lenB = self.molB["len"]

        trace = []

        def trace_from(i, j, k, l):
            for x in self.recursion_cases((i, j, k, l)):
                if guard_case(x[0], (i, j, k, l), self.max_shift):
                    if self.eval_case(x, (i, j, k, l)) == self._M[i, j, k, l]:
                        (io, jo, ko, lo) = x[0]
                        trace.append((io, jo, ko, lo))
                        trace_from(i - io, j - jo, k - ko, l - lo)
                        break

        trace_from(lenA, lenB, lenA, lenB)
        return list(reversed(trace))

    # perform traceback
    # @returns list of 'trace arrows'
    def affine_traceback(self):
        lenA = self.molA["len"]
        lenB = self.molB["len"]

        trace = []

        def shift_by(x, total_shift):
            total_shift[0] += x[0] - x[2]
            total_shift[1] += x[1] - x[3]
            total_shift[2] = abs(total_shift[0]) + abs(total_shift[1])
            return total_shift

        def trace_from(state, idx, total_shift):
            i, j, k, l = idx
            cdef int cidx[4]
            cidx = idx
            if idx == [0, 0, 0, 0] and state == [1, 1, 1, 1]:
                return True
            
            candidates = list()
            for x in self.affine_recursion_cases(state, idx):
                if guard_case(x[1], idx, self.max_shift):
                    if self.affine_eval_case(x, cidx) == self._M[state][idx]:
                        temp_total_shift = total_shift[:]
                        shift_by(x[1], temp_total_shift)
                        shift_by(x[0], temp_total_shift)
                        candidates.append([x,temp_total_shift])

            if candidates:
                sel = argmin([[shift[2],abs(shift[1])] for x,shift in candidates])
                x, temp_total_shift = candidates[sel]
                shift_by(x[1], total_shift)
                (io, jo, ko, lo) = x[1]
                trace.append(x[1])
                return trace_from(x[0], [i - io, j - jo, k - ko, l - lo], total_shift)
            else:    
                return False

        best_score = np.max(
            [ self._M[state][lenA, lenB, lenA, lenB] for state in self.states]
        )
        best_states =[ state for state in self.states 
            if self._M[state][lenA, lenB, lenA, lenB] == best_score ] 

        # select the start state that causes the least number of shifts
        best_state_shifts = [[shift_by(x,[0,0,0])[2]] for x in best_states]
        sel = np.argmin(best_state_shifts)
        best_state = best_states[sel]

        if not trace_from(best_state, [lenA, lenB, lenA, lenB], [0,0,0]):
            print("WARNING: incomplete traceback. Alignment could be garbage.")
        return list(reversed(trace))

    # transfer gap pattern from an alignment string to a sequence string
    @staticmethod
    def _transfer_gaps(alistr, seqstr):
        pos = 0
        res = ""
        for i, c in enumerate(alistr):
            if c == "-":
                res += "-"
            else:
                res += seqstr[pos]
                pos += 1
        return res

    @staticmethod
    def _shift_string(ali, idx):
        length = len(ali[0])

        def shift(i):
            c1 = "X"
            c2 = "X"
            if ali[idx][i] == "-":
                c1 = "-"
            if ali[idx + 2][i] == "-":
                c2 = "-"

            if c1 == c2:
                return "."
            elif c1 == "-":
                return ">"
            elif c2 == "-":
                return "<"

        s = [shift(i) for i in range(length)]
        return "".join(s)

    @staticmethod
    def auto_complete(x, xs):
        xs = list(xs)
        xs.sort()
        for y in xs:
            if y.startswith(x):
                return y
        return x

    # decode trace to alignment strings
    def decode_trace_full(self, trace=None):
        if trace is None:
            trace = self.traceback()

        mols = (self.molA, self.molB, self.molA, self.molB)
        pos = [0] * len(mols)
        alignment = [""] * len(mols)
        for i, y in enumerate(trace):
            for s in range(len(mols)):
                if y[s] == 0:
                    alignment[s] = alignment[s] + "-"
                elif y[s] == 1:
                    alignment[s] = alignment[s] + mols[s]["seq"][pos[s]]
                    pos[s] += 1

        # alignment[0],alignment[1] = highlight_sequence_identity(alignment[0],alignment[1])
        # alignment[2],alignment[3] = highlight_sequence_identity(alignment[2],alignment[3])

        # compute consensus sequences
        cons_seq = [
            consensus_sequence(alignment[2 * i], alignment[2 * i + 1]) for i in range(2)
        ]

        # annotate with structure
        anno_ali = list()
        for alistr, mol in zip(alignment, mols):
            anno_ali.append(self._transfer_gaps(alistr, mol["structure"]))
            anno_ali.append(alistr)

        for i, j in [(4, 6), (0, 2)]:
            if self._is_rna:
                sbpp = consensus_sbpp(
                    alistrA=anno_ali[i],
                    alistrB=anno_ali[j],
                    sbppA=self.molA["sbpp"],
                    sbppB=self.molB["sbpp"],
                )
                structure = mea(sbpp, brackets="[]")[0]
            else:
                structure = consensus_sequence(anno_ali[i], anno_ali[j])
            anno_ali.insert(j + 2, structure)

        shift_strings = list()
        for i in range(2):
            shift_strings.append(self._shift_string(alignment, i))

        alignment = anno_ali

        alignment.insert(len(alignment), cons_seq[1])
        alignment.insert(len(alignment) // 2, cons_seq[0])

        alignment.extend(shift_strings)

        nameA = self._params["nameA"]
        nameB = self._params["nameB"]

        struct_postfix = " ss"
        names = [
            nameA + struct_postfix,
            nameA,
            nameB + struct_postfix,
            nameB,
            "consensus" + struct_postfix,
            "consensus",
            nameA + struct_postfix,
            nameA,
            nameB + struct_postfix,
            nameB,
            "consensus" + struct_postfix,
            "consensus",
            nameA + " shifts",
            nameB + " shifts",
        ]

        return list(zip(names, alignment))

    def decode_trace(self, trace=None):

        alignment = self.decode_trace_full(trace)

        width = max(map(lambda x: len(x[0]), alignment)) + 4

        # add names of single lines
        if not self._params["nodescription"]:
            alignment = [
                ("{:{width}}{}").format(name, alistr, width=width)
                for name, alistr in alignment
            ]
        else:
            alignment = [ alistr for _,alistr in alignment ]

        alignment.append("")

        mode = self.auto_complete(self._params["outmode"], self.outmodes.keys())

        if mode in self.outmodes:
            order = self.outmodes[mode]
        else:
            print(
                "WARNING: unknown output mode. Expect one of "
                + str(list(self.outmodes.keys()))
            )
            order = self.outmodes["sorted"]

        # re-sort
        alignment = [alignment[i] for i in order]

        return alignment

    def eval_affine_trace(self, trace = None):
        cdef int cy[4]
        cdef int cstate[4]

        if trace is None:
            trace = self.traceback()

        def update_state(x,y):
            y = y[:]
            if y[0]==0 and y[1]==0:
                y[0] = x[0]
                y[1] = x[1]
            if y[2]==0 and y[3]==0:
                y[2] = x[2]
                y[3] = x[3]
            return y

        total_score = 0

        beta = self.beta
        gamma = self.gamma
        Delta = self._params["shift_cost"]

        state=[1,1,1,1]
        idx = [0] * 4
        for i, y in enumerate(trace):
            for k in range(4):
                idx[k] += y[k] # <- compute new indices

            i, j, k, l = idx

            mu1 = self.mu1(i, j)
            mu2 = self.mu2(k, l)

            cy = y
            cstate = state

            score = affine_score(cstate, cy, mu1, mu2, beta, gamma, Delta)
            total_score += score

            state = update_state(state,y)
            
            line = " ".join(
                [
                    str(item)
                    for item in [
                        idx,
                        y,
                        score,
                        "-->",
                        total_score,
                        #self._M[state][idx] 
                    ]
                ]
            )
            yield line

    # evaluate trace
    def eval_trace(self, trace = None):
        if self._affine:
            yield from self.eval_affine_trace(trace)
            return

        if trace is None:
            trace = self.traceback()


        idx = [0] * 4
        for i, y in enumerate(trace):
            for k in range(4):
                idx[k] += y[k]
            # lookup case
            for x in self.recursion_cases(idx):
                if x[0] == y:
                    line = " ".join(
                        [
                            str(item)
                            for item in [
                                idx,
                                y,
                                x[1],
                                "-->",
                                self.eval_case(x, idx),
                            ]
                        ]
                    )
                    yield line
                    break


# compute mea structure
def mea(sbpp, gamma=3, *, brackets="()"):
    n = len(sbpp) - 1

    F = np.zeros((n + 1, n + 1), dtype="float")
    T = np.zeros((n + 1, n + 1), dtype="int")

    candidates = [list() for i in range(0, n + 1)]

    for i in reversed(range(1, n + 1)):
        candidates[i].append((i, sbpp[i, i]))
        F[i, i - 1] = 0
        for j in range(i, n + 1):
            for k, C in candidates[j]:
                if F[i, j] < F[i, k - 1] + C:
                    F[i, j] = F[i, k - 1] + C
                    T[i, j] = k

            if i + 3 >= j:
                continue

            C = F[i + 1, j - 1] + 2 * gamma * sbpp[i, j]
            if C > F[i, j]:
                candidates[j].append((i, C))
                F[i, j] = C
                T[i, j] = i

    # trace back

    structure = ["."] * (n + 1)
    stack = list()
    stack.append((1, n))

    while stack:
        (i, j) = stack.pop()
        k = T[i, j]
        if i + 3 >= j or k == 0:
            continue

        if k == j:
            stack.append((i, j - 1))
        elif k == i:
            structure[k] = brackets[0]
            structure[j] = brackets[1]
            stack.append((k + 1, j - 1))
        else:
            stack.append((i, k - 1))
            stack.append((k + 1, j - 1))
            structure[k] = brackets[0]
            structure[j] = brackets[1]

    return ("".join(structure[1:]), F[1, n])


# highlight identical sequence in a pairwise alignment
def highlight_sequence_identity(alistrA, alistrB):
    res = ["", ""]
    for x, y in zip(alistrA.lower(), alistrB.lower()):
        if x == y:
            x = x.upper()
            y = x
        res[0] += x
        res[1] += y
    return res


def consensus_sequence(alistrA, alistrB):
    def c(x, y):
        if x == y:
            return x  # "*"
        else:
            return "."

    return "".join([c(x, y) for x, y in zip(alistrA.upper(), alistrB.upper())])


def parse_dotbracket(dbstr):
    res = [-1] * len(dbstr)
    stack = list()
    for i, sym in enumerate(dbstr):
        if sym == "(":
            stack.append(i)
        elif sym == ")":
            j = stack.pop()
            res[i] = j
            res[j] = i

    return res


# consensus base pair probabilities
def consensus_sbpp(alistrA, sbppA, alistrB, sbppB):
    sbpp = np.zeros((len(alistrA) + 1, len(alistrB) + 1), dtype="float")

    length = [len(sbppA) - 1, len(sbppB) - 1]

    p0 = [1, 1]
    for i0, x0 in enumerate(zip(alistrA, alistrB)):
        p1 = [1, 1]
        for i1, x1 in enumerate(zip(alistrA, alistrB)):
            pr = [0, 0]
            for k, sbppX in [(0, sbppA), (1, sbppB)]:
                if x0[k] == "-" or x1[k] == "-":
                    pr[k] = 0
                else:
                    pr[k] = sbppX[p0[k], p1[k]]

            sbpp[i0 + 1, i1 + 1] = sqrt(pr[0] * pr[1])
            for k in range(2):
                if x1[k] != "-":
                    p1[k] += 1
        for k in range(2):
            if x0[k] != "-":
                p0[k] += 1

    return sbpp


# highlight matched base pairs in a pairwise alignment; for balanced dot bracket strings
def highlight_structure_identity(alistrA, alistrB):

    strA = parse_dotbracket(alistrA)
    strB = parse_dotbracket(alistrB)

    res = ["", ""]
    for i, (x, y) in enumerate(zip(alistrA.lower(), alistrB.lower())):
        if strA[i] >= 0 and strB[i] >= 0 and strA[i] == strB[i]:
            if strA[i] > i:
                x = "["
            else:
                x = "]"
            y = x

        res[0] += x
        res[1] += y

    return res


# highlight matched base pairs in a pairwise alignment; given sbpp matrices
def highlight_structure_similarity(alistrA, alistrB, *, sbppA, sbppB):

    sbpp = consensus_sbpp(alistrA, sbppA, alistrB, sbppB)

    structure = parse_dotbracket(mea(sbpp)[0])

    res = [list(x) for x in [alistrA, alistrB]]
    for i in range(len(alistrA)):
        for j in range(i + 1, len(alistrA)):
            if structure[i] == j:
                res[0][i] = "<"
                res[1][i] = "<"
                res[0][j] = ">"
                res[1][j] = ">"

    return ["".join(x) for x in res]

