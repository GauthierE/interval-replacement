import unittest
import numpy as np

from utils import Representation, Interval


class ComplexTest(unittest.TestCase):

    def setUp(self):
        dimensions = [3, 3]  # 2D grid 3x3
        self.rep = Representation(dimensions)

        self.rep.create_vecs((0, 0), 0)  # (position, vector space's dimension)
        self.rep.create_vecs((0, 1), 1)
        self.rep.create_vecs((0, 2), 1)
        self.rep.create_vecs((1, 0), 2)
        self.rep.create_vecs((1, 1), 2)
        self.rep.create_vecs((1, 2), 2)
        self.rep.create_vecs((2, 0), 2)
        self.rep.create_vecs((2, 1), 2)
        self.rep.create_vecs((2, 2), 2)

        # add linear maps u: x -> y
        self.rep.create_matrix((0, 0), (1, 0), None)  # (x, y, matrix)
        self.rep.create_matrix((0, 0), (0, 1), None)  # for matrices equal to 0 you can directly write None
        self.rep.create_matrix((0, 1), (0, 2), np.array([[0]]))
        self.rep.create_matrix((0, 1), (1, 1), np.array([[0], [1]]))
        self.rep.create_matrix((0, 2), (1, 2), np.array([[0], [1]]))

        self.rep.create_matrix((1, 0), (2, 0), np.array([[0, 0], [0, 0]]))
        self.rep.create_matrix((1, 0), (1, 1), np.array([[0, 0], [0, 0]]))
        self.rep.create_matrix((1, 1), (1, 2), np.array([[1, 0], [0, 0]]))
        self.rep.create_matrix((1, 1), (2, 1), np.array([[1, 0], [0, 1]]))
        self.rep.create_matrix((1, 2), (2, 2), np.array([[1, 0], [0, 1]]))

        self.rep.create_matrix((2, 0), (2, 1), np.array([[1, 0], [0, 1]]))
        self.rep.create_matrix((2, 1), (2, 2), np.array([[1, 0], [0, 0]]))

        self.intervals = self.rep.list_int(conv=False)

    def test_evaluation_from_01_to_21(self):
        assert (self.rep.evaluation((0, 1), (2, 1)) == [[0], [1]]).all()

    def test_evaluation_from_10_to_12(self):
        assert (self.rep.evaluation((1, 0), (1, 2)) == [[0, 0], [0, 0]]).all()

    def test_evaluation_from_11_to_22(self):
        assert (self.rep.evaluation((1, 1), (2, 2)) == [[1, 0], [0, 0]]).all()

    def test_number_of_generated_intervals(self):
        assert len(self.intervals) == 83

    def test_interval_hull_type_3_1(self):
        interval = Interval([(2, 0), (1, 1), (0, 2)], [(2, 2)])
        assert set(self.rep.int_hull(interval)) == {(1, 2), (2, 1), (1, 1), (2, 0), (0, 2), (2, 2)}

    def test_matrix_M_type_3_1(self):
        interval = Interval([(2, 0), (1, 1), (0, 2)], [(2, 2)])
        matrix = self.rep.construct_matrix_M_tot(interval)
        assert (matrix == [[1., 0., -1., -0., 0.],
                           [0., 1., -0., -1., 0.],
                           [1., 0., 0., 0., -0.],
                           [0., 0., 0., 0., -1.],
                           [0., 0., 1., 0., -0.],
                           [0., 0., 0., 0., -1.]]).all()

    def test_matrix_N_type_3_1(self):
        interval = Interval([(2, 0), (1, 1), (0, 2)], [(2, 2)])
        matrix = self.rep.construct_matrix_N_tot(interval)
        assert (matrix == []).all()

    def test_matrix_M_type_2_2(self):
        interval = Interval([(0, 1), (1, 0)], [(2, 1), (1, 2)])
        matrix = self.rep.construct_matrix_M_tot(interval)
        assert (matrix == [[0., -0., -0.],
                           [1., -0., -0.]]).all()

    def test_matrix_N_type_2_2(self):
        interval = Interval([(0, 1), (1, 0)], [(2, 1), (1, 2)])
        matrix = self.rep.construct_matrix_N_tot(interval)
        assert (matrix == [[1.,  0.],
                           [0.,  1.],
                           [-1., -0.],
                           [-0., -0.]]).all()

    def test_interval_rank_type_3_1(self):
        interval = Interval([(2, 0), (1, 1), (0, 2)], [(2, 2)])
        assert self.rep.int_rank(interval) == 0

    def test_interval_replacement_type_3_1(self):
        interval = Interval([(2, 0), (1, 1), (0, 2)], [(2, 2)])
        assert self.rep.int_replacement(interval) == 0

    def test_interval_replacement_all(self):
        interval_replacement = {}
        for itv in self.intervals:
            repl = self.rep.int_replacement(itv)
            if repl != 0:
                interval_replacement[(tuple(sorted(itv.src)), tuple(sorted(itv.snk)))] = repl

        assert interval_replacement == {(((1, 0),), ((1, 0),)): 2,
                                        (((0, 2),), ((2, 2),)): 1,
                                        (((0, 1), (2, 0)), ((2, 1),)): 1,
                                        (((1, 1), (2, 0)), ((2, 2),)): 1}
