class SWAlign(object):
    """SW alignment for two sequence."""
    equal_val = 3
    unequ_val = -3
    gap_value = -2

    def __init__(self, query, target, match_score=3, mis_score=2,
                 gap_score=-1):
        """Init class."""
        self.query = query
        self.target = target
        self.query_len = len(query)
        self.target_len = len(target)
        self.match_score = match_score
        self.mis_score = mis_score
        self.gap_score = gap_score
        self.alignment()

    def alignment(self):
        """Get maxtrix for alignment."""
        self.align_score = None
        self.align_pos = []
        self.maxtrix = []
        for i in range(0, self.query_len + 1):
            tmp = []
            for j in range(0, self.target_len + 1):
                tmp.append(0)
            self.maxtrix.append(tmp)
        self.max_info = []
        max_val = 0
        for i in range(1, self.query_len + 1):
            for j in range(1, self.target_len + 1):
                tmpvalue = 0
                if self.query[i-1] == self.target[j-1]:
                    tmpvalue = self.equal_val
                else:
                    tmpvalue = self.unequ_val
                self.maxtrix[i][j] = max(
                    self.maxtrix[i-1][j-1] + tmpvalue,
                    self.maxtrix[i-1][j] + self.gap_value,
                    self.maxtrix[i][j-1] + self.gap_value, 0)
                if self.maxtrix[i][j] > max_val:
                    max_val = self.maxtrix[i][j]
                    self.max_info = []
                    self.max_info.append((i, j))
                elif self.maxtrix[i][j] == max_val:
                    self.max_info.append((i, j))
        align_info = {}
        for tmp in self.max_info:
            x, y = tmp
            pos = []
            while x >= 1 and y >= 1:
                maxvalue = max(self.maxtrix[x-1][y-1], self.maxtrix[x-1][y],
                               self.maxtrix[x][y-1])
                if maxvalue == 0:
                    break
                else:
                    if maxvalue == self.maxtrix[x-1][y-1]:
                        pos.append((x, y))
                        x -= 1
                        y -= 1
                    elif maxvalue == self.maxtrix[x-1][y]:
                        pos.append((x, 0))
                        x -= 1
                    else:
                        pos.append((0, y))
                        y -= 1
            pos.reverse()
            score = 0
            for xval, yval in pos:
                if xval != 0 and yval != 0:
                    if self.query[xval-1] == self.target[yval-1]:
                        score += self.match_score
                    else:
                        score += self.mis_score
                else:
                    score += self.gap_score
            align_info.setdefault(score, []).append(pos)
        maxscore = max(align_info.keys())
        self.align_score = maxscore
        self.align_pos = align_info[maxscore][0]

    def show_alignment(self, pos):
        """Show alignment sequence."""
        def get_seq(seq, tmppos):
            if tmppos == 0:
                return '-'
            else:
                try:
                    return seq[tmppos-1]
                except IndexError:
                    return '-'
        xstart, ystart = 0, 0
        xseq, yseq = '', ''
        for (x, y) in pos:
            xseq += get_seq(self.query, x)
            yseq += get_seq(self.target, y)
            if xstart == 0:
                xstart = x
            if ystart == 0:
                ystart = y
        return xstart, ystart, xseq, yseq

    def disapply_align(self, disapply_len=20):
        """Disapply alignment."""
        pre_len = len(str(max(x for x in self.align_pos[-1])))
        header_len = 8 + pre_len
        header_format = '{:-%ss}' % (header_len)
        print(header_format)
        print(pre_len)
        for i in range(0, len(self.align_pos), disapply_len):
            tmppos = self.align_pos[i:i+disapply_len]
            print(tmppos)
            xstart, ystart, xseq, yseq = self.show_alignment(tmppos)
            line1_header = 'Query: ' + str(xstart)
            print(line1_header)
            line1 = header_format.format(line1_header) + xseq
            line2 = ' ' * header_len
            for x, y in zip(xseq, yseq):
                showtag = ' '
                if x != '-' and y != '-':
                    if x == y:
                        showtag = '|'
                line2 += showtag
            line3_header = 'Sbjct: ' + str(ystart)
            line3 = header_format.format(line3_header) + yseq
            yield [line1, line2, line3]


def test():
    """This is a test case."""
    seq1 = 'TGTTACGG'
    seq2 = 'GGTTGACTA'
    sw = SWAlign(seq1, seq2)
    print(sw.maxtrix)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')


if __name__ == '__main__':
    test()
