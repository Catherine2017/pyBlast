"""
This is Needleman–Wunsch algorithm.

See more on https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm.
"""

from align import Align


class NWalign(Align):
    """NW alignment for two sequence."""

    def __init__(self, query, target, subs_matrix='../matrices/NUC.4.4',
                 gap_penatly=(-2,)):
        """Init class."""
        super().__init__(query, target, subs_matrix, gap_penatly)
        if len(gap_penatly) == 1:
            self.do_alignment_linear(gap_penatly[0])
        elif len(gap_penatly) == 2:
            self.do_alignment_affine(*gap_penatly)

    def do_alignment_linear(self, gap_value):
        """Do NW alignment by linear gap penalty."""
        # 初始化矩阵
        self.maxtrix = []
        for i in range(0, self.target_len + 1):
            tmp = []
            for j in range(0, self.query_len + 1):
                val = None
                if i == 0:
                    val = j * gap_value
                elif j == 0:
                    val = i * gap_value
                tmp.append(val)
            self.maxtrix.append(tmp)
        # 计算得分
        self.traback = {}
        for i in range(1, self.target_len + 1):
            for j in range(1, self.query_len + 1):
                tmpvalue = self.subs_matrix[self.target[i-1]][self.query[j-1]]
                self.maxtrix[i][j] = max(
                    self.maxtrix[i-1][j-1] + tmpvalue,
                    self.maxtrix[i-1][j] + gap_value,
                    self.maxtrix[i][j-1] + gap_value
                )
                if self.maxtrix[i][j] == self.maxtrix[i-1][j-1] + tmpvalue:
                    self.traback[(i, j)] = (i-1, j-1)
                elif self.maxtrix[i][j] == self.maxtrix[i-1][j] + gap_value:
                    self.traback[(i, j)] = (i-1, j)
                else:
                    self.traback[(i, j)] = (i, j-1)
        # 回溯发现最优比对
        self.align_score = self.maxtrix[self.target_len][self.query_len]
        x, y = self.target_len, self.query_len
        align_query, align_target = [], []
        while x > 0 and y > 0:
            x1, y1 = self.traback[(x, y)]
            if x1 == x - 1 and y1 == y - 1:
                align_target.append(self.get_seq(self.target, x))
                align_query.append(self.get_seq(self.query, y))
                self.align_pos.append([x, y])
            elif x1 == x - 1:
                align_target.append(self.get_seq(self.target, x))
                align_query.append('-')
                self.align_pos.append([x, 0])
            else:
                align_target.append('-')
                align_query.append(self.get_seq(self.query, y))
                self.align_pos.append([0, y])
            x = x1
            y = y1
        self.align_pos.reverse()
        self.align_query = ''.join(align_query[::-1])
        self.align_target = ''.join(align_target[::-1])

    def do_alignment_affine(self, gap_open, gap_exten):
        """Do NW alignment by affine gap penalty."""
        pass


def test():
    """This is a test case."""
    print('-' * 24)
    seq1 = 'TACGGGCCCGCTAC'
    seq2 = 'TAGCCCTATCGGTCA'
    sw = NWalign(seq1, seq2)
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
    sw = NWalign(seq1, seq2, '../matrices/BLOSUM62')
    print(sw.align_pos)
    print(sw.align_query)
    print(sw.align_target)
    print(sw.align_score)
    for lines in sw.disapply_align():
        print('\n'.join(lines))
        print('\n')


if __name__ == '__main__':
    test()
