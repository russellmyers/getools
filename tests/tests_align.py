class AlignTests(unittest.TestCase):

    def test_simple(self):

        seq1 = Sequence('panamabananas')
        seq2 = Sequence('anamaqqbananas')
        self.assertTrue(seq1.is_valid())
        al = GlobalAligner()
        score = al.align(seq1, seq2)
        self.assertEqual(score, 9)
        res = al.backtrack()
        self.assertCountEqual(res, ['panama--bananas', '-anamaqqbananas'])

        seq1 = DNASequence('AGCTTAGCT')
        seq2 = DNASequence('ACCGAGCT')
        self.assertTrue(seq1.is_valid())
        al = GlobalAligner(scoring_matrix=ScoringMatrix(DNASequence.SIMPLE_SCORING_MATRIX, DNASequence.DNA_ALPHABET))
        score = al.align(seq1, seq2)
        self.assertEqual(score, 2)
        res = al.backtrack()
        self.assertCountEqual(res,['AGCTTAGCT', 'A-CCGAGCT'] )

