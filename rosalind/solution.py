import numpy as np
import collections
import math
import re
from collections import defaultdict, Counter, OrderedDict
import sys
sys.setrecursionlimit(10000000)


nucleotides = "ACGT"
comp_nucleotides = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def compliment(dna):
    return ''.join([comp_nucleotides[n] for n in dna])


def comp_rev(dna):
    return compliment(dna)[::-1]


def mutate(s, d):
    if d == 0:
        return [s]
    if len(s) == 0:
        return [""]
    m = mutate(s[1:], d)
    q = mutate(s[1:], d - 1)
    t = []
    for i in m:
        t.append(s[0] + i)
    for i in q:
        for j in nucleotides:
            t.append(j + i)
    return t


def hamming_dist(a, b):
    assert len(a) == len(b), "{}, {}".format(a, b)
    return sum([0 if a[i] == b[i] else 1 for i in range(len(a))])


def k_mers(k, s):
    for i in range(len(s) - k + 1):
        yield s[i: i+k]


def count_occurrence(s, p, d=0):
    ls, lp = len(s), len(p)
    cnt = 0
    for i in range(ls - lp + 1):
        if hamming_dist(p, s[i: i+lp]) <= d:
            cnt += 1
    return cnt


def p1(inpfile="", oupfile=""):
    print(open(inpfile, "r").read().replace('T', 'U'))


def p2(inpfile="", oupfile=""):
    DNA, k, _ = open(inpfile, "r").read().split('\n')
    k = int(k)
    d = {}
    for i in range(len(DNA)):
        if DNA[i: i+k] not in d:
            d[DNA[i: i+k]] = 0
        d[DNA[i: i+k]] += 1
    ml = 0
    for i in d:
        if d[i] > ml:
            ml = d[i]

    for i in d:
        if d[i] == ml:
            print(i, end=" ")


def p3(inpfile="", oupfile=""):
    inputs = open(inpfile, "r").read().split('\n')
    s = inputs[0]
    k, d = map(int, inputs[1].strip('\n').split(' '))
    cnt = defaultdict(int)
    for i in range(len(s) - k + 1):
        t = set(mutate(s[i: i+k], d))
        for dna in t:
            cnt[dna] += 1

    kmer = list(cnt.keys())
    ml = max([cnt[i] + cnt[comp_rev(i)] for i in kmer])

    for i in cnt:
        if cnt[i] + cnt[comp_rev(i)] == ml:
            print(i, end=" ")
            # print(i, count_occurrence(s, i, d), count_occurrence(s, comp_rev(i), d), ml)


def p3_norev(inpfile="", oupfile=""):
    inputs = open(inpfile, "r").read().split('\n')
    s = inputs[0]
    k, d = map(int, inputs[1].strip('\n').split(' '))
    cnt = defaultdict(int)
    for i in range(len(s) - k + 1):
        t = set(mutate(s[i: i+k], d))
        for dna in t:
            cnt[dna] += 1

    kmer = list(cnt.keys())
    ml = max([cnt[i] for i in kmer])

    for i in cnt:
        if cnt[i] == ml:
            print(i, end=" ")
            # print(i, count_occurrence(s, i, d), count_occurrence(s, comp_rev(i), d), ml)


def p5(inpfile="", oupfile=""):
    print(hamming_dist(*tuple(open(inpfile, "r").read().strip('\n').split('\n'))))


class InpScanner(object):
    def __init__(self, inp):
        self.inps = [t for t in inp.split()]
        self.lines = [t for t in inp.split('\n')]
        self.cnt = 0
        self.tot = len(self.inps)

    def next(self):
        self.cnt += 1
        return self.inps[self.cnt - 1]

    def next_int(self):
        return int(self.next())

    def next_float(self):
        return float(self.next())

    def eof(self):
        return self.cnt == self.tot

    def iter(self):
        while not self.eof():
            yield self.next()


class RosalindSolver(object):
    """
    by default ouput are put to stdout
    """
    def __init__(self, inpid='sample'):
        inpfile = "rosalind_{}.txt".format(inpid)
        self.inp = InpScanner(open(inpfile, "r").read())

    def solve(self):
        raise NotImplemented


class P4(RosalindSolver):
    def solve(self):
        dna = self.inp.next()
        ld = len(dna)
        skew = np.cumsum([0] + [1 if dna[i] is 'G' else (-1 if dna[i] is 'C' else 0) for i in range(ld)])
        print(*np.where(skew == skew.min())[0].tolist())


def dist(patt, dna):
    if isinstance(dna, list):
        return sum([dist(patt, d) for d in dna])
    else:
        return min([hamming_dist(patt, k_mer) for k_mer in k_mers(len(patt), dna)])


class P6(RosalindSolver):
    """
    I implement the solution of complexity 4^k * n * t * k
    """

    @staticmethod
    def enumerate(k):
        """
        can use itertools.product
        :param k:
        :return:
        """
        def decode(i):
            return ''.join([nucleotides[(i >> j) % 4] for j in range(0, k * 2, 2)])

        for i in range(4 ** k):
            yield decode(i)

    def solve(self):
        k = self.inp.next_int()
        inps = [i for i in self.inp.iter()]

        ans = ""
        min_dist = 23333333
        for k_mer in self.enumerate(k):
            temp_dist = dist(k_mer, inps)
            if temp_dist < min_dist:
                min_dist = temp_dist
                ans = k_mer

        print(ans)


def prob_by_prof(prof, k_mer):
    prob = 1
    for i in range(len(k_mer)):
        prob *= prof[nucleotides.index(k_mer[i])][i]
    return prob


class P7(RosalindSolver):
    def solve(self):
        dna = self.inp.next()
        k = self.inp.next_int()
        prof = [[self.inp.next_float() for i in range(k)] for i in range(4)]
        ans = ""
        max_prob = 0

        for k_mer in k_mers(k, dna):
            temp_prob = prob_by_prof(prof, k_mer)
            if temp_prob > max_prob:
                max_prob = temp_prob
                ans = k_mer
        print(ans)


def profile(motifs, smoothing=True):
    t = len(motifs)
    sf = 1 if smoothing else 0
    c2p = lambda cnt: [(cnt[n] + sf) / (t + sf * 4) for n in nucleotides]
    return list(zip(*[c2p(Counter(s)) for s in zip(*motifs)]))


def consensus(motifs):
    return ''.join([Counter(s).most_common(1)[0][0] for s in zip(*motifs)])


def score(motifs):
    return dist(consensus(motifs), motifs)
    # def count(s):
    #     return sum(sorted([s.count(n) for n in nucleotides])[:3])
    # return sum([count(s) for s in zip(*motifs)])


class P8(RosalindSolver):
    """
    In P8, I implement the GibbsSampler, but it have to run 2 mins
    """
    def solve(self):
        k, t, N = self.inp.next_int(), self.inp.next_int(), self.inp.next_int()
        dna = [self.inp.next() for i in range(t)]
        ld = len(dna[0])
        for i in range(t):
            assert ld == len(dna[i])

        niter = 1
        best_motifs = None
        min_score = 233333

        def offset2motifs(offset):
            return [dna[i][offset[i]: offset[i] + k] for i in range(t)]

        def profile_sample(prof, text):
            def normalize(arr):
                tot = sum(arr)
                return [p / tot for p in arr]

            p = [prob_by_prof(prof, k_mer) for k_mer in k_mers(k, text)]
            i = np.random.choice(list(range(ld - k + 1)), p=normalize(p))
            return text[i: i+k]

        for iter in range(niter):
            motifs = offset2motifs(np.random.randint(0, ld - k + 1, size=(t,)))
            temp_score = score(motifs)
            if temp_score < min_score:
                min_score = temp_score
                best_motifs = motifs
            for j in range(N):
                i = np.random.randint(0, t)
                prof = profile(motifs[:i] + motifs[i+1:])
                motifs[i] = profile_sample(prof, dna[i])
                temp_score = score(motifs)
                if temp_score < min_score:
                    min_score = temp_score
                    best_motifs = motifs

        for motif in best_motifs:
            print(motif)


class P9(RosalindSolver):
    def solve(self):
        k = self.inp.next_int()
        dna = self.inp.next()
        for t in k_mers(k, dna):
            print(t)


class Arc(object):
    def __init__(self, fr, to, label, next):
        self.fr, self.to, self.label = fr, to, label
        self.next = next
        self.visited = False


class AdjacentLink(object):
    def __init__(self):
        self.h = defaultdict(lambda: -1)
        self.cin = defaultdict(lambda: 0)
        self.cou = defaultdict(lambda: 0)
        self.arc = []

    def add(self, fr, to, label):
        self.cin[to] += 1
        self.cou[fr] += 1
        self.arc.append(Arc(fr, to, label, self.h[fr]))
        self.h[fr] = len(self.arc) - 1

    def out_arcs(self, cur):
        i = self.h[cur]
        while i != -1:
            yield self.arc[i]
            i = self.arc[i].next


class P10(RosalindSolver):
    def solve(self):
        k = self.inp.next_int()
        texts = [t for t in self.inp.iter()]
        adj = AdjacentLink()

        for t in texts:
            fr, to, label = t[:-1], t[1:], t
            adj.add(fr, to, label)

        def eulerian_path(cur):
            paths = []
            for a in adj.out_arcs(cur):
                if not a.visited:
                    a.visited = True
                    paths.append(eulerian_path(a.to))
            return sum(paths, [cur])

        paths = None
        for t in texts:
            fr, to, label = t[:-1], t[1:], t
            if adj.cou[fr] - adj.cin[fr]  == 1:
                paths = eulerian_path(fr)
                break

        print(''.join([paths[0]] + [p[-1] for p in paths[1:]]))


class P11(RosalindSolver):
    def solve(self):
        texts = [t for t in self.inp.iter()]
        adj = AdjacentLink()
        for t in texts:
            fr, to, label = t[:-1], t[1:], t
            adj.add(fr, to, label)

        for fr in set([t[:-1] for t in texts]):
            print("{} -> {}".format(fr, ','.join([a.to for a in adj.out_arcs(fr)])))


class P12(RosalindSolver):
    def solve(self):
        texts = [t for t in self.inp.iter()]
        adj = AdjacentLink()

        for t in texts:
            fr, to, label = t[:-1], t[1:], t
            adj.add(fr, to, label)

        ver = set([t[:-1] for t in texts])
        non_branching_ver = set([v for v in ver if adj.cin[v] == adj.cou[v] == 1])
        # print(ver, non_branching_ver)
        # tails = []
        # for u in non_branching_ver:
        #     for v in adj.out_arcs(u):
        #         tails.append(v)
        # heads = non_branching_ver.difference(set(tails))

        def path_walking(u):
            if u not in non_branching_ver:
                return [u]
            for a in adj.out_arcs(u):
                return [u] + path_walking(a.to)

        for u in ver.difference(non_branching_ver):
            for a in adj.out_arcs(u):
                path = [u] + path_walking(a.to)
                print(path[0] + ''.join([t[-1] for t in path[1:]]), end=' ')


codon_table = {
            "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V", "UUC": "F", "CUC": "L", "AUC": "I",
            "GUC": "V", "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V", "UUG": "L", "CUG": "L",
            "AUG": "M", "GUG": "V", "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A", "UCC": "S",
            "CCC": "P", "ACC": "T", "GCC": "A", "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
            "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A", "UAU": "Y", "CAU": "H", "AAU": "N",
            "GAU": "D", "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D", "UAA": "Stop", "CAA": "Q",
            "AAA": "K", "GAA": "E", "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E", "UGU": "C",
            "CGU": "R", "AGU": "S", "GGU": "G", "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
            "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G", "UGG": "W", "CGG": "R", "AGG": "R",
            "GGG": "G"
        }


def translate(rna):
    l = len(rna)
    return ''.join([codon_table[rna[i: i+3]] for i in range(0, l, 3)])


def transcript(dna):
    return dna.replace('T', 'U')


class P13(RosalindSolver):
    def solve(self):
        dna = self.inp.next()
        peptide = self.inp.next()

        k = len(peptide) * 3

        for kmer in k_mers(k, dna):
            if translate(transcript(kmer)) == peptide or translate(transcript(comp_rev(kmer))) == peptide:
                print(kmer)


class P14(RosalindSolver):
    def solve(self):
        s1 = list(map(float, self.inp.lines[0].split()))
        s2 = list(map(float, self.inp.lines[1].split()))
        # print(len(s1), len(s2), len([i - j for j in s2 for i in s1]))
        minkowski_diff = Counter([round(np.abs(i - j), 5) for j in s2 for i in s1])
        ans = minkowski_diff.most_common(1)[0]
        print("{}\n{}\n".format(ans[1], ans[0]))


amino_acid = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115,
               'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


def mass_of_peptide(pep):
    return sum([amino_acid[i] for i in pep])


def linear_spectrum(peptide):
    prefix_mass = [0]
    lp = len(peptide)
    for i in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid[i])
    res = [0]
    for i in range(0, lp):
        for j in range(i + 1, lp + 1):
            res.append(prefix_mass[j] - prefix_mass[i])
    return sorted(res)


def cyclic_spectrum(peptide):
    prefix_mass = [0]
    lp = len(peptide)
    peptide = peptide + peptide
    for i in peptide:
        prefix_mass.append(prefix_mass[-1] + amino_acid[i])
    res = [0, prefix_mass[lp]]
    for i in range(0, lp):
        for j in range(1, lp):
            res.append(prefix_mass[i + j] - prefix_mass[i])

    return sorted(res)


class P15(RosalindSolver):
    def solve(self):
        print(*linear_spectrum(self.inp.next()), end=' ')


class P16(RosalindSolver):
    def solve(self):
        spectrum = sorted([int(v) for v in self.inp.iter()])

        parent_mass = max(spectrum)
        peptides = set([""])
        ans = []
        while len(peptides) != 0:
            peptides = set([t + n for n in amino_acid for t in peptides])
            removed = []
            # print(peptides)
            for pep in peptides:
                if mass_of_peptide(pep) == parent_mass:
                    pep_spectrum = cyclic_spectrum(pep)
                    # print(len(pep_spectrum), len(spectrum))
                    if pep_spectrum == spectrum:
                        ans.append(pep)
                    removed.append(pep)
                else:
                    t = linear_spectrum(pep)
                    for i in t:
                        if i not in spectrum:
                            removed.append(pep)
                            break
            peptides = peptides.difference(set(removed))

        def print_ans(pep):
            return '-'.join([str(mass_of_peptide(t)) for t in pep])

        print(' '.join(set([print_ans(pep) for pep in ans])))


class P17(RosalindSolver):
    def solve(self):
        a = self.inp.next()
        b = self.inp.next()
        la, lb = len(a), len(b)
        f = [[0 for j in range(lb + 1)] for i in range(la + 1)]
        p = [[0 for j in range(lb + 1)] for i in range(la + 1)]
        for i in range(la):
            for j in range(lb):
                if f[i - 1][j] > f[i][j - 1]:
                    p[i][j] = -1
                    f[i][j] = f[i - 1][j]
                else:
                    p[i][j] = 1
                    f[i][j] = f[i][j - 1]

                if a[i] == b[j] and f[i - 1][j - 1] + 1 > f[i][j]:
                    f[i][j] = f[i - 1][j - 1] + 1
                    p[i][j] = 0

        def print_path(i, j):
            if i == -1 or j == -1:
                return []
            if p[i][j] == 0:
                return print_path(i - 1, j - 1) + [a[i]]
            elif p[i][j] == -1:
                return print_path(i - 1, j)
            else:
                return print_path(i, j - 1)

        # print(f[la - 1][lb - 1])
        print(''.join(print_path(la - 1, lb - 1)))


class P18(RosalindSolver):
    def solve(self):
        _ = self.inp.next()
        s, t = "", ""
        for seg in self.inp.iter():
            if seg.startswith('>'):
                break
            s += seg
        for seg in self.inp.iter():
            t += seg

        ls, lt = len(s), len(t)
        f = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        p = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        for i in range(ls):
            for j in range(lt):
                if f[i - 1][j] < f[i][j - 1]:
                    p[i][j] = 1
                    f[i][j] = f[i - 1][j] + 1
                else:
                    p[i][j] = -1
                    f[i][j] = f[i][j - 1] + 1

                if s[i] == t[j]:
                    if f[i - 1][j - 1] < f[i][j]:
                        f[i][j] = f[i - 1][j - 1]
                        p[i][j] = 0
                elif f[i - 1][j - 1] + 1 < f[i][j]:
                        f[i][j] = f[i - 1][j - 1] + 1
                        p[i][j] = 0

        def print_path(i, j):
            if i == -1 and j == -1:
                return [], []
            if i == -1:
                return ['-'] * j, [c for c in t[:j]]
            if j == -1:
                return [c for c in s[:i]], ['-'] * i

            if p[i][j] == 0:
                a, b = print_path(i - 1, j - 1)
                return a + [s[i]], b + [t[j]]
            elif p[i][j] == -1:
                a, b = print_path(i, j - 1)
                return a + ['-'], b + [t[j]]
            elif p[i][j] == 1:
                a, b = print_path(i - 1, j)
                return a + [s[i]], b + ['-']


        # print(f[la - 1][lb - 1])
        a, b = print_path(ls - 1, lt - 1)
        print("{}\n{}\n{}\n".format(f[ls - 1][lt - 1], ''.join(a), ''.join(b)))


class P19(RosalindSolver):
    def solve(self):
        score = {
            'A': {'A': 4, 'C': 0, 'D': -2, 'E': -1, 'F': -2, 'G': 0, 'H': -2, 'I': -1, 'K': -1, 'L': -1,
                  'M': -1, 'N': -2, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 0, 'V': 0, 'W': -3, 'Y': -2},
            'C': {'A': 0, 'C': 9, 'D': -3, 'E': -4, 'F': -2, 'G': -3, 'H': -3, 'I': -1, 'K': -3, 'L': -1,
                  'M': -1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -1, 'T': -1, 'V': -1, 'W': -2, 'Y': -2},
            'D': {'A': -2, 'C': -3, 'D': 6, 'E': 2, 'F': -3, 'G': -1, 'H': -1, 'I': -3, 'K': -1, 'L': -4,
                  'M': -3, 'N': 1, 'P': -1, 'Q': 0, 'R': -2, 'S': 0, 'T': -1, 'V': -3, 'W': -4, 'Y': -3},
            'E': {'A': -1, 'C': -4, 'D': 2, 'E': 5, 'F': -3, 'G': -2, 'H': 0, 'I': -3, 'K': 1, 'L': -3,
                  'M': -2, 'N': 0, 'P': -1, 'Q': 2, 'R': 0, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
            'F': {'A': -2, 'C': -2, 'D': -3, 'E': -3, 'F': 6, 'G': -3, 'H': -1, 'I': 0, 'K': -3, 'L': 0,
                  'M': 0, 'N': -3, 'P': -4, 'Q': -3, 'R': -3, 'S': -2, 'T': -2, 'V': -1, 'W': 1, 'Y': 3},
            'G': {'A': 0, 'C': -3, 'D': -1, 'E': -2, 'F': -3, 'G': 6, 'H': -2, 'I': -4, 'K': -2, 'L': -4,
                  'M': -3, 'N': 0, 'P': -2, 'Q': -2, 'R': -2, 'S': 0, 'T': -2, 'V': -3, 'W': -2, 'Y': -3},
            'H': {'A': -2, 'C': -3, 'D': -1, 'E': 0, 'F': -1, 'G': -2, 'H': 8, 'I': -3, 'K': -1, 'L': -3,
                  'M': -2, 'N': 1, 'P': -2, 'Q': 0, 'R': 0, 'S': -1, 'T': -2, 'V': -3, 'W': -2, 'Y': 2},
            'I': {'A': -1, 'C': -1, 'D': -3, 'E': -3, 'F': 0, 'G': -4, 'H': -3, 'I': 4, 'K': -3, 'L': 2,
                  'M': 1, 'N': -3, 'P': -3, 'Q': -3, 'R': -3, 'S': -2, 'T': -1, 'V': 3, 'W': -3, 'Y': -1},
            'K': {'A': -1, 'C': -3, 'D': -1, 'E': 1, 'F': -3, 'G': -2, 'H': -1, 'I': -3, 'K': 5, 'L': -2,
                  'M': -1, 'N': 0, 'P': -1, 'Q': 1, 'R': 2, 'S': 0, 'T': -1, 'V': -2, 'W': -3, 'Y': -2},
            'L': {'A': -1, 'C': -1, 'D': -4, 'E': -3, 'F': 0, 'G': -4, 'H': -3, 'I': 2, 'K': -2, 'L': 4,
                  'M': 2, 'N': -3, 'P': -3, 'Q': -2, 'R': -2, 'S': -2, 'T': -1, 'V': 1, 'W': -2, 'Y': -1},
            'M': {'A': -1, 'C': -1, 'D': -3, 'E': -2, 'F': 0, 'G': -3, 'H': -2, 'I': 1, 'K': -1, 'L': 2,
                  'M': 5, 'N': -2, 'P': -2, 'Q': 0, 'R': -1, 'S': -1, 'T': -1, 'V': 1, 'W': -1, 'Y': -1},
            'N': {'A': -2, 'C': -3, 'D': 1, 'E': 0, 'F': -3, 'G': 0, 'H': 1, 'I': -3, 'K': 0, 'L': -3,
                  'M': -2, 'N': 6, 'P': -2, 'Q': 0, 'R': 0, 'S': 1, 'T': 0, 'V': -3, 'W': -4, 'Y': -2},
            'P': {'A': -1, 'C': -3, 'D': -1, 'E': -1, 'F': -4, 'G': -2, 'H': -2, 'I': -3, 'K': -1, 'L': -3,
                  'M': -2, 'N': -2, 'P': 7, 'Q': -1, 'R': -2, 'S': -1, 'T': -1, 'V': -2, 'W': -4, 'Y': -3},
            'Q': {'A': -1, 'C': -3, 'D': 0, 'E': 2, 'F': -3, 'G': -2, 'H': 0, 'I': -3, 'K': 1, 'L': -2,
                  'M': 0, 'N': 0, 'P': -1, 'Q': 5, 'R': 1, 'S': 0, 'T': -1, 'V': -2, 'W': -2, 'Y': -1},
            'R': {'A': -1, 'C': -3, 'D': -2, 'E': 0, 'F': -3, 'G': -2, 'H': 0, 'I': -3, 'K': 2, 'L': -2,
                  'M': -1, 'N': 0, 'P': -2, 'Q': 1, 'R': 5, 'S': -1, 'T': -1, 'V': -3, 'W': -3, 'Y': -2},
            'S': {'A': 1, 'C': -1, 'D': 0, 'E': 0, 'F': -2, 'G': 0, 'H': -1, 'I': -2, 'K': 0, 'L': -2,
                  'M': -1, 'N': 1, 'P': -1, 'Q': 0, 'R': -1, 'S': 4, 'T': 1, 'V': -2, 'W': -3, 'Y': -2},
            'T': {'A': 0, 'C': -1, 'D': -1, 'E': -1, 'F': -2, 'G': -2, 'H': -2, 'I': -1, 'K': -1, 'L': -1,
                  'M': -1, 'N': 0, 'P': -1, 'Q': -1, 'R': -1, 'S': 1, 'T': 5, 'V': 0, 'W': -2, 'Y': -2},
            'V': {'A': 0, 'C': -1, 'D': -3, 'E': -2, 'F': -1, 'G': -3, 'H': -3, 'I': 3, 'K': -2, 'L': 1,
                  'M': 1, 'N': -3, 'P': -2, 'Q': -2, 'R': -3, 'S': -2, 'T': 0, 'V': 4, 'W': -3, 'Y': -1},
            'W': {'A': -3, 'C': -2, 'D': -4, 'E': -3, 'F': 1, 'G': -2, 'H': -2, 'I': -3, 'K': -3, 'L': -2,
                  'M': -1, 'N': -4, 'P': -4, 'Q': -2, 'R': -3, 'S': -3, 'T': -2, 'V': -3, 'W': 11, 'Y': 2},
            'Y': {'A': -2, 'C': -2, 'D': -3, 'E': -2, 'F': 3, 'G': -3, 'H': 2, 'I': -1, 'K': -2, 'L': -1,
                  'M': -1, 'N': -2, 'P': -3, 'Q': -1, 'R': -2, 'S': -2, 'T': -2, 'V': -1, 'W': 2, 'Y': 7}}

        sigma, epsilon = -11, -1

        s, t = self.inp.next(), self.inp.next()
        ls, lt = len(s), len(t)
        f = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        u = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        l = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        pf = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        pu = [[0 for j in range(lt + 1)] for i in range(ls + 1)]
        pl = [[0 for j in range(lt + 1)] for i in range(ls + 1)]

        for i in range(ls):
            for j in range(lt):
                if u[i][j - 1] + epsilon < f[i][j - 1] + sigma:
                    u[i][j] = f[i][j - 1] + sigma
                    pu[i][j] = 1
                else:
                    u[i][j] = u[i][j - 1] + epsilon
                    pu[i][j] = 0

                if l[i - 1][j] + epsilon < f[i - 1][j] + sigma:
                    l[i][j] = f[i - 1][j] + sigma
                    pl[i][j] = 1
                else:
                    l[i][j] = l[i - 1][j] + epsilon
                    pl[i][j] = 0

                if u[i][j] < l[i][j]:
                    f[i][j] = l[i][j]
                    pf[i][j] = -1
                else:
                    f[i][j] = u[i][j]
                    pf[i][j] = 1

                if f[i][j] < f[i - 1][j - 1] + score[s[i]][t[j]]:
                    f[i][j] = f[i - 1][j - 1] + score[s[i]][t[j]]
                    pf[i][j] = 0

        def print_path(i, j, p=0):
            if i == -1 and j == -1:
                return [], []
            if i == -1:
                return ['-'] * j, [c for c in t[:j]]
            if j == -1:
                return [c for c in s[:i]], ['-'] * i

            if p == -1 and pl[i + 1][j] == 0:
                a, b = print_path(i - 1, j, -1)
                return a + [s[i]], b + ['-']
            if p == 1 and pu[i][j + 1] == 0:
                a, b = print_path(i, j - 1, 1)
                return a + ['-'], b + [t[j]]
            if pf[i][j] == 0:
                a, b = print_path(i - 1, j - 1)
                return a + [s[i]], b + [t[j]]
            elif pf[i][j] == 1:
                a, b = print_path(i, j - 1, 1)
                return a + ['-'], b + [t[j]]
            elif pf[i][j] == -1:
                a, b = print_path(i - 1, j, -1)
                return a + [s[i]], b + ['-']


        # print(f[la - 1][lb - 1])
        a, b = print_path(ls - 1, lt - 1)
        print("{}\n{}\n{}\n".format(f[ls - 1][lt - 1], ''.join(a), ''.join(b)))


class P20(RosalindSolver):
    def solve(self):
        nums = [t for t in self.inp.iter()]
        nums[0] = nums[0][1:]
        nums[-1] = nums[-1][:-1]
        nums = list(map(int, nums))
        ans = 0
        for i in range(len(nums) - 1):
            if nums[i + 1] - nums[i] != 1:
                ans += 1
        print(ans + ans % 2)


class P21(RosalindSolver):
    def solve(self):
        chrs = re.findall('\(.*?\)', self.inp.lines[0], re.S)
        edges = []
        for chr in chrs:
            nodes = list(map(int, re.findall('[+-]\d+', chr)))
            last = nodes[-1] * 2 if nodes[-1] > 0 else (-nodes[-1] * 2) - 1
            for n in nodes:
                if n < 0:
                    edges.append((last, -n * 2))
                    last = (-n * 2) - 1
                else:
                    edges.append((last, n * 2 - 1))
                    last = n * 2
        print(', '.join(["({}, {})".format(*e) for e in edges]))


class P22(RosalindSolver):
    def solve(self):
        nodes = re.findall('\d+', self.inp.lines[0])
        cnt = len(nodes)
        edges = [(int(nodes[i]), int(nodes[i + 1])) for i in range(0, cnt, 2)]
        index = {int(nodes[i]): i // 2 for i in range(0, cnt, 2)}
        # print(edges, index)
        _visited = [0 for e in edges]

        def visited(e):
            return _visited[index[e[0]]] == 1

        def set_visited(e):
            _visited[index[e[0]]] = 1

        def visit(e):
            if visited(e):
                return []
            set_visited(e)
            if e[1] % 2 == 0:
                cur = -e[1] // 2
                e = edges[index[(-cur * 2) - 1]]
            else:
                cur = (e[1] + 1) // 2
                e = edges[index[cur * 2]]
            return [cur] + visit(e)

        ans = []
        for e in edges:
            if visited(e):
                continue
            else:
                ans.append(visit(e))
        print(''.join(['(' + ' '.join(['{:+d}'.format(d) for d in cycle]) + ')' for cycle in ans]))


class P23(RosalindSolver):
    def solve(self):
        inp = ''.join(self.inp.inps)
        inps = re.findall('[ATGC]+', inp)
        n = len(inps)
        m = len(inps[0])
        ans = [[round(hamming_dist(a, b) / m, 5) for a in inps] for b in inps]
        for row in ans:
            for i in row:
                print(i, end=' ')
            print()


class P24(RosalindSolver):
    def solve(self):
        lines = [line for line in self.inp.lines if len(line) > 0]
        tot_line = len(lines)

        x, y = None, None

        def split(t):
            stack = []
            res = []
            cur = []
            for c in t + ',':
                if c == ',':
                    if len(stack) == 0:
                        res.append(''.join(cur))
                        cur = []
                        continue
                cur.append(c)
                if c == '(':
                    stack.append(c)
                if c == ')':
                    del stack[-1]
            return res

        def solve(desc, dep=0):
            if '(' not in desc:
                return dep
            sub, this = re.match('\((.*)\)(.*)', desc).groups()
            # import ipdb
            # ipdb.set_trace()
            # print(sub[:10], this, x, y)
            subdescs = split(sub)
            dx, dy = None, None
            if this == x:
                dx = dep
            if this == y:
                dy = dep
            for subdesc in subdescs:
                if x in subdesc and y in subdesc:
                    return solve(subdesc, dep+1)
                if x in subdesc:
                    dx = solve(subdesc, dep+1)
                if y in subdesc:
                    dy = solve(subdesc, dep+1)
            if x in desc and y in desc:
                return dx + dy - dep * 2
            if dx is not None:
                return dx
            return dy

        for i in range(0, tot_line, 2):
            desc = lines[i]
            x, y = lines[i+1].split()
            print(solve(desc[:-1]), end=' ')




class P25(RosalindSolver):
    def solve(self):
        n = self.inp.next_int()
        dis = np.array([[self.inp.next_float() for i in range(n)] for j in range(n)], dtype=float)

        ans = []

        def neighbor_join(D, n, step, remaining):
            if n == 2:
                ans.append((remaining[0], remaining[1], D[0, 1]))
                print(ans[-1])
                return
            _D = (n - 2) * D - D.sum(axis=0, keepdims=True) - D.sum(axis=1, keepdims=True)
            np.fill_diagonal(_D, np.inf)
            (i, j) = np.unravel_index(np.argmin(_D), _D.shape)
            if i > j:
                i, j = j, i
            np.fill_diagonal(_D, 0)
            print(_D)
            totdis = D.sum(axis=0)
            delta = (totdis[i] - totdis[j]) / (n - 2)
            d = D[i, j]
            limlen_i = (d + delta) / 2
            limlen_j = (d - delta) / 2
            newD = np.zeros([n + 1, n + 1])
            remaining.append(step)
            ans.append((step, remaining[i], limlen_i))
            print(ans[-1])
            ans.append((step, remaining[j], limlen_j))
            print(ans[-1])
            newD[:-1, :-1] = D
            for k in range(n):
                newD[n, k] = newD[k, n] = 0.5 * (D[k, i] + D[k, j] - D[i, j])
            newD = np.delete(newD, [i, j], axis=0)
            newD = np.delete(newD, [i, j], axis=1)
            del remaining[i]
            del remaining[j - 1]
            neighbor_join(newD, n - 1, step + 1, remaining)

        neighbor_join(dis, n, n, list(range(n)))
        for t in ans:
            print("{}->{}:{}".format(*t))
            print("{}->{}:{}".format(t[1], t[0], t[2]))


class P26(RosalindSolver):
    def solve(self):
        k, m = self.inp.next_int(), self.inp.next_int()
        inp_tmp = [float(t) for t in self.inp.iter()]
        p = [inp_tmp[i: i+m] for i in range(0, len(inp_tmp), m)]
        n = len(p)

        def dist(a, b):
            return sum([(a[i] - b[i]) ** 2 for i in range(m)])

        center = p[:k]

        def closest(p):
            x = 0
            d = dist(p, center[0])
            for i in range(1, k):
                _d = dist(p, center[i])
                if _d < d:
                    d = _d
                    x = i
            return x

        def average(cluster):
            a = [0.0 for i in range(m)]
            for c in cluster:
                a = [a[i] + c[i] for i in range(m)]
            nc = len(cluster)
            a = [a[i] / nc for i in range(m)]
            return a

        def iterate():
            assign = [closest(t) for t in p]
            for i in range(k):
                center[i] = average([p[j] for j in range(n) if assign[j] == i])

        for i in range(1110):
            iterate()

        for t in center:
            print(" ".join(["{}".format(round(t[i], 3)) for i in range(m)]))


class P27(RosalindSolver):
    def solve(self):
        n = self.inp.next_int()
        dis = [[self.inp.next_float() for i in range(n)] for j in range(n)]
        clusters = [[i] for i in range(n)]

        def dist(a, b):
            return sum([dis[i][j] for i in a for j in b]) / len(a) / len(b)

        for _ in range(n - 1):
            d = float('inf')
            a, b = None, None
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    newd = dist(clusters[i], clusters[j])
                    if newd < d:
                        d = newd
                        a, b = i, j
            clusters[a] += clusters[b]
            del clusters[b]
            print(" ".join([str(i + 1) for i in sorted(clusters[a])]))


class P28(RosalindSolver):
    def solve(self):
        text = self.inp.next()
        n = len(text)
        sorted_suffix = sorted(text[i:] for i in range(n))
        ans = [n - len(suf) for suf in sorted_suffix]
        print(', '.join([str(a) for a in ans]))


class P29(RosalindSolver):
    def solve(self):
        text = self.inp.next()
        n = len(text)
        sorted_suffix = sorted(text[i:] for i in range(n))
        ans = [n - len(suf) for suf in sorted_suffix]
        print(''.join([text[a - 1] for a in ans]))


class P30(RosalindSolver):
    def solve(self):
        text = self.inp.next()
        n = len(text)
        ans = ["" for i in range(n)]
        for i in range(n):
            ans = [text[j] + ans[j] for j in range(n)]
            ans = sorted(ans)
        print(ans[0][1:] + ans[0][0])


if __name__ == "__main__":
    # p1("rosalind_rna.txt")
    # p2("rosalind_ba1b.txt")
    # p3("rosalind_ba1j.txt")
    # p3_norev("rosalind_ba1i.txt")
    # P4('ba1f').solve()
    # p5("rosalind_ba1g.txt")
    # P6('ba2b').solve()
    # P7('ba2c').solve()
    # P8('ba2g').solve()
    # P9('ba3a').solve()
    # P10('ba3h').solve()
    # P11('ba3e').solve()
    # P12('ba3k').solve()
    # P13('ba4b').solve()
    # P14('conv').solve()
    # P15('ba4j').solve()
    # P16('ba4e').solve()
    # P17('ba5c').solve()
    # P18('edta').solve()
    # P19('ba5j').solve()
    # P20('ba6b').solve()
    # P21('ba6h').solve()
    # P22('ba6i').solve()
    # P23('pdst').solve()
    P24('nwck').solve()
    # P25('ba7e').solve()
    # P25().solve()
    # P26('ba8c').solve()
    # P27('ba8e').solve()
    # P28('ba9g').solve()
    # P29('ba9i').solve()
    # P30('ba9j').solve()

