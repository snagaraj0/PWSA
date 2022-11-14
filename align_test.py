import unittest

from align import *

TEST_INPUT_FILE="test_example.input"

class TestAlignmentClasses(unittest.TestCase):

    def test_match_matrix(self):
        """
        Tests match matrix object
        """
        match_matrix = MatchMatrix()
        match_matrix.set_score("A", "C", 5)
        self.assertEqual(match_matrix.get_score("A", "C"), 5)

    def test_score_matrix_score(self):
        """
        Tests score matrix object score set + get methods
        """
        score_matrix = ScoreMatrix("M", 2, 2)
        score_matrix.set_score(1, 1, 5.3)
        self.assertEqual(score_matrix.get_score(1, 1), 5.3)

    def test_score_matrix_pointers(self):
        """
        Tests score matrix object pointer set + get methods
        """
        score_matrix = ScoreMatrix("M", 2, 2)
        score_matrix.set_pointers(1, 1, score_matrix.score_matrix[0][0])
        score_matrix.set_pointers(1, 1, score_matrix.score_matrix[0][1])
        score_matrix.set_pointers(1, 1, score_matrix.score_matrix[1][0])
        self.assertEqual(score_matrix.get_pointers(1, 1), [("M", 0, 0), ("M", 0, 1), ("M", 1, 0)])
    def test_param_loading(self):
        """
        Tests AlignmentParameters "load_params_from_file()" function
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)
        self.assertEqual(align_params.seq_a, "AATGC")
        self.assertEqual(align_params.seq_b, "AGGC")
        self.assertTrue(align_params.global_alignment)
        self.assertEqual(align_params.dx, 0.1)
        self.assertEqual(align_params.ex, 0.5)
        self.assertEqual(align_params.dy, 0.6)
        self.assertEqual(align_params.ey, 0.3)
        self.assertEqual(align_params.alphabet_a, "ATGC")
        self.assertEqual(align_params.alphabet_b, "ATGCX")
        self.assertEqual(align_params.len_alphabet_a, 4)
        self.assertEqual(align_params.len_alphabet_b, 5)

        # test that the match matrix is set up correctly
        #  if this fails, make sure you are loading the asymmetric matrix properly!
        match_mat = align_params.match_matrix
        self.assertEqual(match_mat.get_score("A", "X"), 0.3)
        self.assertEqual(match_mat.get_score("C", "G"), -0.3)
        self.assertEqual(match_mat.get_score("G", "C"), 0)


    def test_update_ix(self):
        """
        Test AlignmentAlgorithm's update Ix
        """

        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.ix_matrix.set_score(2, 2, 2.5)

        # run the method!
        align.update_ix(3, 2)

        score = align.ix_matrix.get_score(3,2)
        self.assertEqual(score, 2)

        # Check pointers
        self.assertEqual(align.ix_matrix.get_pointers(3,2), [("M", 2, 2), ("Ix", 2, 2)])


    def test_update_m(self):
        """
        Test AlignmentAlgorithm's update M
        """
        # Create match matrix
        align_params = AlignmentParameters()
        match = MatchMatrix()
        match.set_score("C", "A", 0)
        match.set_score("A", "T", -1)
        match.set_score("A", "C", -1)
        match.set_score("G", "C", -1)

        align_params.match_matrix = match
        align_params.seq_a = "AGAC"
        align_params.seq_b = "CCTA"
        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 1, 2.2)
        align.ix_matrix.set_score(2, 1, 1.8)
        align.iy_matrix.set_score(2, 1, 1.8)

        align.update_m(3, 2)
        score = align.m_matrix.get_score(3, 2)
        self.assertEqual(round(score, 1), 1.2)
    
    def test_update_iy(self):
        """
        Test AlignmentAlgorithm's update Iy
        """
        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dx = 1
        align_params.ex = 0.5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 3)
        align.iy_matrix.set_score(2, 2, 2.5)

        # run the method!
        align.update_iy(2, 3)

        score = align.iy_matrix.get_score(2, 3)
        self.assertEqual(score, 2)


    def test_traceback_start(self):
        """
        Tests that the traceback finds the correct start
        Should test local and global alignment!
        """

        # Global traceback test
        align_params = AlignmentParameters()
        align_params.global_alignment = True

        align = Align("", "")
        align.align_params = align_params
        align.m_matrix = ScoreMatrix("M", 3, 3)
        for row in range(align.m_matrix.nrow):
            for col in range(align.m_matrix.ncol):
                align.m_matrix.score_matrix[row][col].score = 2
        align.m_matrix.score_matrix[2][1].score = 6
        self.assertEqual(align.find_traceback_start(), (6, {(2,1)}))


        # Local traceback test
        align.m_matrix.score_matrix[1][2].score = 8
        self.assertEqual(align.find_traceback_start(), (8, {(1,2)}))




if __name__=='__main__':
    unittest.main()
