"""
Align for two sequence.

See more on https://en.wikipedia.org/wiki/Sequence_alignment.
"""
import os
import re


class Align(object):
    """Show alignment for two sequence."""

    def __init__(self, query, target, subs_matrix, gap_penatly):
        """Init class."""
        self.load_subsmatrix(subs_matrix)
        checkset = set(self.subs_matrix.keys())
        self.query = self.check_sequence(query, checkset)
        self.target = self.check_sequence(target, checkset)
        self.query_len = len(query)
        self.target_len = len(target)
        self.align_score = 0
        self.align_pos = []
        self.align_query = ''
        self.align_target = ''
        self.gap_penatly = gap_penatly

    def load_subsmatrix(self, matrixfile):
        """Load substitution matrix file."""
        self.subs_matrix = {}
        with open(matrixfile, 'rt') as rd:
            cols = []
            for line in rd:
                line = line.rstrip(os.linesep)
                if line.startswith(('#', '\n')):
                    continue
                item = line.strip().split()
                if re.search(r'^\s+', line):
                    cols = item
                else:
                    for i in range(1, len(item)):
                        self.subs_matrix.setdefault(item[0], {}).setdefault(
                            cols[i-1], int(item[i]))

    def check_sequence(self, sequence, checkset):
        """Check sequence type."""
        sequence = re.sub(r'\s+', '', sequence).upper()
        aset = set(sequence)
        if not aset.issubset(set(checkset)):
            raise TypeError("Sequence %s is not subset of %s!" % (
                aset, checkset))
        return sequence

    def show_alignment(self, pos):
        """Show alignment sequence."""
        xstart, ystart = 0, 0
        for (x, y) in pos:
            if xstart == 0:
                xstart = x
            if ystart == 0:
                ystart = y
        return xstart, ystart

    def disapply_align(self, disapply_len=40):
        """Disapply alignment."""
        pre_len = len(str(max(x for x in self.align_pos[-1])))
        header_len = 8 + pre_len
        header_format = '{:<%ss}' % (header_len)
        for i in range(0, len(self.align_pos), disapply_len):
            tmppos = self.align_pos[i:i+disapply_len]
            yseq = self.align_query[i:i+disapply_len]
            xseq = self.align_target[i:i+disapply_len]
            xstart, ystart = self.show_alignment(tmppos)
            line1_header = 'Query: ' + str(ystart)
            line1 = header_format.format(line1_header) + yseq
            line2 = ' ' * header_len
            for x, y in zip(xseq, yseq):
                showtag = ' '
                if x != '-' and y != '-':
                    if x == y:
                        showtag = '|'
                line2 += showtag
            line3_header = 'Sbjct: ' + str(xstart)
            line3 = header_format.format(line3_header) + xseq
            yield [line1, line2, line3]

    def get_seq(self, seq, tmppos):
        """Get pos unit."""
        if tmppos >= 1:
            return seq[tmppos-1]
        else:
            return '-'

    def get_score(self, pos):
        """Calculate score for alignment."""
        score = 0
        gap_count = 0
        align_query, align_target = '', ''
        for xval, yval in pos:
            xunit = self.get_seq(self.target, xval)
            yunit = self.get_seq(self.query, yval)
            align_query += yunit
            align_target += xunit
            if xval != 0 and yval != 0:
                score += self.subs_matrix[xunit][yunit]
                gap_count = 0
            else:
                gap_count += 1
                if gap_count > 1 and len(self.gap_penatly) > 1:
                        score += self.gap_penatly[1]
                else:
                    score += self.gap_penatly[0] + self.gap_penatly[1]
        return score, align_query, align_target
