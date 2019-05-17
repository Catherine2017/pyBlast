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
        super().__init__(query, target, subs_matrix)
        if len(gap_penatly) == 1:
            self.do_alignment_linear(gap_penatly[0])
        elif len(gap_penatly) == 2:
            self.do_alignment_affine(*gap_penatly)

    def do_alignment_linear(self, gap_value):
        """Do NW alignment by linear gap penalty."""
        # 初始化矩阵
        self.maxtrix = []
        

    def do_alignment_affine(gap_open, gap_exten):
        """Do NW alignment by affine gap penalty."""
        pass
