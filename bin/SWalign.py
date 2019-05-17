"""
This is Smith–Waterman algorithm.

See more on https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm.
"""

from align import Align


class SWAlign(Align):
    """SW alignment for two sequence."""

    def __init__(self, query, target, subs_matrix='../matrices/NUC.4.4',
                 gap_penatly=(-2,)):
        """Init class."""
        super().__init__(query, target, subs_matrix, gap_penatly)
        if len(gap_penatly) == 1:
            self.do_alignment_linear(gap_penatly[0])
        elif len(gap_penatly) == 2:
            self.do_alignment_affine(*gap_penatly)

    def do_alignment_linear(self, gap_value):
        """Do SW alignment by linear gap penalty."""
        self.maxtrix = []
        # 初始化矩阵
        for i in range(0, self.target_len + 1):
            tmp = []
            for j in range(0, self.query_len + 1):
                val = None
                if i == 0 or j == 0:
                    val = 0
                tmp.append(val)
            self.maxtrix.append(tmp)
        max_info = []
        self.traback = {}
        max_val = 0
        # 获得矩阵和最大值对应的坐标
        for i in range(1, self.target_len + 1):
            for j in range(1, self.query_len + 1):
                tmpvalue = self.subs_matrix[self.target[i-1]][self.query[j-1]]
                maxtmp = max(
                    self.maxtrix[i-1][j-1] + tmpvalue,
                    self.maxtrix[i-1][j] + gap_value,
                    self.maxtrix[i][j-1] + gap_value)
                if maxtmp == self.maxtrix[i-1][j-1] + tmpvalue:
                    self.traback[(i, j)] = (i-1, j-1)
                elif maxtmp == self.maxtrix[i-1][j] + gap_value:
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
                x = x1
                y = y1
                if self.maxtrix[x][y] == 0:
                    break
            pos.reverse()
            """
            if pos[0][0] == 0 or pos[0][1] == 0:
                del pos[0]"""
            # 计算比对分值
            score, align_query, align_target = self.get_score(pos)
            align_info.setdefault(score, []).append({
                'align_pos': pos, 'align_query': align_query,
                'align_target': align_target})
        # 找到得分最高的比对
        maxscore = max(align_info.keys())
        self.align_score = maxscore
        self.align_pos = align_info[maxscore][0]['align_pos']
        self.align_query = align_info[maxscore][0]['align_query']
        self.align_target = align_info[maxscore][0]['align_target']

    def do_alignment_affine(gap_open, gap_exten):
        """Do SW alignment by affine gap penalty."""
        pass


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
    seq1 = 'MFLSILVALCLWLHLALGVRGAPCEAVRIPMCRHMPWNITRMPNHLHHSTQENAILAIEQ'
    seq2 = 'MFLSILVALCLWLHLALGVRGAPCEAVRICRHMPWNITRMPNHLNAILAIEQ'
    sw = SWAlign(seq1, seq2, '../matrices/BLOSUM62')
    print(sw.align_pos)
    print(sw.align_query)
    print(sw.align_target)
    print(sw.align_score)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')
    seq1 = 'CGCTGCGGGAGCGGAGCGGGTCGGTGCGGCCGG'
    seq2 = 'CGCTGCGGGAGCGGCTGCCGGGGFFGGTGCGGCCGG'
    sw = SWAlign(seq1, seq2)


if __name__ == '__main__':
    test()
