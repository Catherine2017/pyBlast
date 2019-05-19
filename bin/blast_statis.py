"""Blast statistics for alignment."""
from math import exp, factorial, log


class BlastStat(object):
    """Blast statistics for alignment."""

    def calculate_entropy(self, sequence, logval=2):
        """Shannon Entropy Calculator."""
        tot, adict = 0, {}
        for i in sequence:
            tot += 1
            adict[i] = adict.setdefault(i, 0) + 1
        H = 0
        for i, val in adict.items():
            per = val/tot
            H += per * log(per)
        H = -H/log(logval)
        return H

    def calcuate_lambda(self, match, mismatch, pn=0.25):
        """Estimating lambda."""
        expect_score = match * pn + mismatch * (1-pn)
        lamb, high, low = 1, 2, 0
        while high - low > 0.001:
            sum = pn * pn * exp(lamb * match) * 4\
                + pn * pn * exp(lamb * mismatch) * 12
            if sum > 1:
                high = lamb
                lamb = (lamb + low)/2
            else:
                low = lamb
                lamb = (lamb + high)/2
        targetID = pn * pn * exp(lamb * match) * 4
        H = lamb * match * targetID + lamb * mismatch * (1 - targetID)
        return expect_score, lamb, H, targetID


class OneHSP(object):
    # http://etutorials.org/Misc/blast/Part+III+Practice/Chapter+7.+A+BLAST+Statistics+Tutorial/7.1+Basic+BLAST+Statistics/

    def effe_query_length(self, query_len, hsp_len, gapk):
        m_prime = query_len - hsp_len
        if m_prime < 1/gapk:
            return 1/gapk
        else:
            return m_prime

    def effe_db_length(self, dblen, hsp_len, num_seqs, gapk):
        n_prime = dblen - (num_seqs*hsp_len)
        if n_prime < 1/gapk:
            return 1/gapk
        else:
            return n_prime

    @staticmethod
    def n2b(n):
        return n/log(2)

    @staticmethod
    def b2n(n):
        return n*log(2)

    def rawscore2bitscore(self, raw_score, gapk, lamb):
        return self.n2b(lamb*raw_score - log(gapk))

    def rawscore2expect(self, raw_score, k, lamb, m, n):
        # m: effe_query_length
        # n: effe_db_length
        return k*m*n*exp(-1*lamb*raw_score)

    def bitscore2expect(self, bitscore, m, n):
        # m: effe_query_length
        # n: effe_db_length
        return m*n*2**(-1*bitscore)


class MoreHSP(object):
    """Sum Statistics Are Pair-Wise in Their Focus."""

    def sumScore(self, raw_scores, k, lamb, m, n, g):
        # m # effective length of query sequence
        # n # effective length of sbjct sequence
        # g # gap_size;for NCBI-BLAST this value is 50
        r = len(raw_scores)
        if r <= 1:
            raise ValueError("Single score!")
        total_raw_scores = sum(raw_scores)
        n_score = lamb * total_raw_scores
        return n_score - log(k*m*n) - (r - 1) * (
            log(k) + 2*log(g) - log(factorial(r)))
