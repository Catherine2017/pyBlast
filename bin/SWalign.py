"""
This is Smith–Waterman algorithm.
See more on https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm.
"""

from align import Align


class SWAlign(Align):
    """SW alignment for two sequence."""

    equal_val = 3
    unequ_val = -3
    gap_value = -2

    def __init__(self, query, target):
        """Init class."""
        self.query = query
        self.target = target
        self.query_len = len(query)
        self.target_len = len(target)
        self.alignment()

    def alignment(self):
        """Get maxtrix for alignment."""
        self.align_score = None
        self.align_pos = []
        self.maxtrix = []
        # 初始化矩阵
        for i in range(0, self.target_len + 1):
            tmp = []
            for j in range(0, self.query_len + 1):
                tmp.append(0)
            self.maxtrix.append(tmp)
        max_info = []
        self.traback = {}
        max_val = 0
        # 获得矩阵和最大值对应的坐标
        for i in range(1, self.target_len + 1):
            for j in range(1, self.query_len + 1):
                tmpvalue = 0
                if self.target[i-1] == self.query[j-1]:
                    tmpvalue = self.equal_val
                else:
                    tmpvalue = self.unequ_val
                maxtmp = max(
                    self.maxtrix[i-1][j-1] + tmpvalue,
                    self.maxtrix[i-1][j] + self.gap_value,
                    self.maxtrix[i][j-1] + self.gap_value)
                if maxtmp == self.maxtrix[i-1][j-1] + tmpvalue:
                    self.traback[(i, j)] = (i-1, j-1)
                elif maxtmp == self.maxtrix[i-1][j] + self.gap_value:
                    self.traback[(i, j)] = (i-1, j)
                else:
                    self.traback[(i, j)] = (i, j-1)
                self.maxtrix[i][j] = max(maxtmp, 0)
                if self.maxtrix[i][j] > max_val:
                    max_val = self.maxtrix[i][j]
                    max_info = []
                    max_info.append((i, j))
                elif self.maxtrix[i][j] == max_val:
                    max_info.append((i, j))
        align_info = {}
        print(self.maxtrix)
        # 回溯发现比对路径
        for tmp in max_info:
            x, y = tmp
            pos = []
            while x >= 1 and y >= 1:
                x1, y1 = self.traback[(x, y)]
                print(self.maxtrix[x][y])
                if x1 == x - 1 and y1 == y - 1:
                    pos.append((x, y))
                elif x1 == x - 1:
                    pos.append((x, 0))
                else:
                    pos.append((0, y))
<<<<<<< HEAD
=======
                x = x1
                y = y1
>>>>>>> 1e4db8ca882df9a71ee7afe258882788521520c0
                print(self.maxtrix[x][y])
                if self.maxtrix[x][y] == 0:
                    break
            pos.reverse()
            """
            if pos[0][0] == 0 or pos[0][1] == 0:
                del pos[0]"""
            # 计算比对分值
            score = 0
            align_query, align_target = '', ''
            for xval, yval in pos:
                xunit = self.get_seq(self.target, xval)
                yunit = self.get_seq(self.query, yval)
                align_query += yunit
                align_target += xunit
                if xval != 0 and yval != 0:
                    if xunit == yunit:
                        score += self.equal_val
                    else:
                        score += self.unequ_val
                else:
                    score += self.gap_value
            align_info.setdefault(score, []).append({
                'align_pos': pos, 'align_query': align_query,
                'align_target': align_target})
        # 找到得分最高的比对
        maxscore = max(align_info.keys())
        self.align_score = maxscore
        self.align_pos = align_info[maxscore][0]['align_pos']
        self.align_query = align_info[maxscore][0]['align_query']
        self.align_target = align_info[maxscore][0]['align_target']

    @staticmethod
    def get_seq(seq, tmppos):
        """Get pos unit."""
        if tmppos >= 1:
            return seq[tmppos-1]
        else:
            return '-'

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


def test():
    """This is a test case."""
    print('-' * 24)
    seq1 = 'TGTTACGG'
    seq2 = 'GGTTGACTA'
    sw = SWAlign(seq1, seq2)
    print(sw.align_pos)
    print(sw.align_query)
    print(sw.align_target)
    print(sw.align_score)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')
    print('-' * 24)
    seq1 = 'TACGGGCCCGCTAC'
    seq2 = 'TAGCCCTATCGGTCA'
    sw = SWAlign(seq1, seq2)
    print(sw.align_pos)
    print(sw.align_query)
    print(sw.align_target)
    print(sw.align_score)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')
    print('-' * 24)
    seq1 = 'CGCTGCGGGAGCGGAGCGGGTCGGTGCGGCCGG'
    seq2 = 'CGCTGCGGGAGCGGCTGCCGGGGTCGGTGCGGCCGG'
    sw = SWAlign(seq1, seq2)
    print(sw.align_pos)
    print(sw.align_query)
    print(sw.align_target)
    print(sw.align_score)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')


if __name__ == '__main__':
    test()
