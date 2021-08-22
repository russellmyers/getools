import unittest
from getools.utils.string_utils import BWT, pattern_match_naive, pattern_match_zbox, CountingSort, RadixSort,\
     CountingSortExp, CountingSortExpSpecialList, SuffixArray, gen_rand_string, SuffixArrayRankList, RadixSortSpecialList
from getools.utils.file_utils import file_to_list, read_txt_file, read_fasta
from getools.popdist import PopDist, PopDistGen
from getools.cross import ChromosomeTemplate, GenomeTemplate, Organism
import random
import numpy as np
import math

class BWTTests(unittest.TestCase):
  def test_panama(self):
      s = 'panamabananas'
      bw = BWT(s)
      self.assertEqual(bw.bw, 'smnpbnnaaaaa$a')
      orig = bw.reconstruct_naive(also_return_suffix_array=True)
      self.assertEqual(len(orig), 2)
      self.assertEqual(orig[0],s)
      self.assertEqual(len(orig[1]), len(s)+1)
      self.assertEqual(orig[1][4], 7)


  def test_panama_match(self):
        s = 'panamabananas'
        bw = BWT(s)
        res = bw.pattern_match('ana')
        self.assertCountEqual(res, [1, 7, 9])

  def test_panama_from_bw(self):
      s = 'panamabananas'
      bw = BWT(s)
      bw_2 = BWT.from_bw(bw.bw, bw.sa)
      res = bw.pattern_match('ana')
      self.assertCountEqual(res, [1, 7, 9])

  def test_no_suffix_ar(self):
      bw = BWT.from_bw('enwvpeoseu$llt')
      res1 = bw.pattern_match('one')
      self.assertCountEqual(res1, ['bwtind:7']) #No real positions as no suffix array
      res2  = bw.pattern_match('one')
      self.assertCountEqual(res2, ['bwtind:7']) # Still no real positions as force_reconstruct hasnt been used
      res3 = bw.pattern_match('plus', force_reconstruct=True)
      self.assertCountEqual(res3, [6])
      res4  = bw.pattern_match('one') #Now returns real positions as force reconstruct has been used
      self.assertCountEqual(res4, [10])

  def test_ecoli_from_file(self):
      sequences = read_fasta('../test_data/large/GCF_000005845.2_ASM584v2_genomic.fna')
      s = sequences[0][1][:500000]
      #t = s[100000:101000]  # 'GATGATGAT'
      t = 'GATTACA'
      print(len(s), len(t), len(s) * len(t), len(s) + len(t))
      matches = pattern_match_naive(s, t)
      expected_matches = [23254, 80864, 155458, 181444, 262735, 313586, 316375, 322262, 357439, 386155, 392966, 393371, 479364, 481239, 498552]
      self.assertCountEqual(matches, expected_matches)
      matches = pattern_match_zbox(s, t)
      self.assertCountEqual(matches, expected_matches)
      bw_s = read_txt_file('../test_data/large/ecoli_500000.bw')
      bw_sa = file_to_list('../test_data/large/ecoli_500000.sa')
      bw_sa = [int(x) for x in bw_sa]
      bw = BWT.from_bw(bw_s, sa=bw_sa)
      # bw = BWT(s)
      matches = bw.pattern_match(t)
      self.assertCountEqual(matches, expected_matches)

  def test_ecoli_construct_bwt(self):
      sequences = read_fasta('../test_data/large/GCF_000005845.2_ASM584v2_genomic.fna')
      s = sequences[0][1][:50000]
      #t = s[100000:101000]  # 'GATGATGAT'
      t = 'GATTACA'
      print(len(s), len(t), len(s) * len(t), len(s) + len(t))
      matches = pattern_match_naive(s, t)
      expected_matches = [23254] #[42, 485, 701, 996, 1007, 1178, 1379, 1745, 2094, 2170, 2387, 2609, 3336, 3394, 3880, 4054, 4285, 4371, 4398, 4486, 4896, 4945, 4978, 5143, 5149, 5792, 5926, 5951, 6526, 6836, 7317, 7380, 8368, 8389, 9670, 9934, 9999, 10047, 10469, 10681, 11204, 11447, 11778, 11789, 11852, 11910, 12212, 12365, 12377, 12503, 12546, 14179, 14190, 14512, 14690, 14928, 15124, 15430, 15772, 15873, 15916, 15958, 16008, 16036, 16098, 16350, 16515, 16757, 17094, 17143, 17199, 17300, 17419, 17431, 17504, 17565, 17802, 17847, 18273, 18378, 18396, 18558, 18833, 19298, 19357, 19547, 19584, 19655, 19675, 19752, 19868, 19920, 19963, 19991, 20153, 20340, 20612, 20659, 20812, 21061...
      self.assertCountEqual(matches, expected_matches)
      bw = BWT(s, verbose=1)
      matches = bw.pattern_match(t)
      self.assertCountEqual(matches, expected_matches)



  def test_zbox(self):
      sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAAAAAAATTAGGCGATGCGGAGC'
      res, Z = pattern_match_zbox(sample_s, 'AAG', also_return_zbox=True)
      exp_Z = [-1, 1, 0, 0, 2, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 2, 2, 2, 2, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0]
      self.assertEqual(res, [2])
      self.assertEqual(Z, exp_Z)
      _, Z = pattern_match_zbox("aaaaaa", None, also_return_zbox=True)  # X 5 4 3 2 1
      self.assertEqual(Z, [-1, 5, 4, 3, 2, 1])
      _, Z = pattern_match_zbox("aabaacd", None, also_return_zbox=True)  # X 1 0 2 1 0 0
      self.assertEqual(Z, [-1, 1, 0, 2, 1, 0, 0])
      _, Z = pattern_match_zbox("abababab", None, also_return_zbox=True)  # X 0 6 0 4 0 2 0
      self.assertEqual(Z, [-1, 0, 6, 0, 4, 0, 2, 0])
      res, Z = pattern_match_zbox("abababab", "bab", also_return_zbox=True)  # X 0 6 0 4 0 2 0
      self.assertCountEqual(res, [1,3,5])
      self.assertEqual(Z, [-1, 0, 1, 0, 0, 3, 0, 3, 0, 3, 0, 1])

      _, Z = pattern_match_zbox("ACACACGTACACG", None, also_return_zbox=True)  # X 0 4 0 2 0 0 0 4 0 2 0 0
      self.assertEqual(Z,[-1, 0, 4, 0, 2, 0, 0, 0, 4, 0, 2, 0, 0] )


class CountingSortTests(unittest.TestCase):

    def test_simple(self):
        cs = CountingSort('banana')
        res = cs.sort()
        self.assertEqual(res,'aaabnn')
        cs = CountingSort([5, 12, 0, -3, 15, 12,  2])
        res = cs.sort()
        self.assertEqual(res, [-3, 0, 2, 5, 12, 12,15])
        cs = CountingSort(['b', 'a', 'n', 'a', 'n', 'a'])
        res = cs.sort()
        self.assertEqual(res, ['a', 'a', 'a', 'b', 'n', 'n'])

    def test_edge(self):
        cs = CountingSort('')
        res = cs.sort()
        self.assertEqual(res,'')
        cs = CountingSort([])
        res = cs.sort()
        self.assertEqual(res,[])



    def test_misc(self):
        cs = CountingSort('banana')
        res = cs.sort()
        self.assertEqual(res, 'aaabnn')
        cs = CountingSort([5, 12, 0, -3, 15, 2])
        res = cs.sort()
        self.assertEqual(res,[-3, 0, 2, 5, 12, 15])
        cs = CountingSort(['b', 'a', 'n', 'a', 'n', 'a'])
        res = cs.sort()
        self.assertEqual(res, ['a', 'a', 'a', 'b', 'n', 'n'])
        cs = CountingSort('5032671932')
        res = cs.sort()
        self.assertEqual(res, '0122335679')


class CountingSortExpTests(unittest.TestCase):

    def test_string(self):
        with self.assertRaises(Exception) as context:
            cs = CountingSortExp('banana')
            res = cs.sort()
        self.assertTrue("Error - Must be list. Type provided is: <class 'str'>" in str(context.exception))


    def test_list_of_strings(self):
        cs = CountingSortExp(['abc', 'a', 'bcae', 'bcad', 'ab'], exp=3)
        res = cs.sort()
        self.assertEqual(res,['abc', 'a', 'ab', 'bcad', 'bcae'])

        cs = CountingSortExp(['abc', 'a', 'bcae', 'bcad', 'ab'], exp=0)
        res = cs.sort()
        self.assertEqual(res,['abc', 'a', 'ab', 'bcae', 'bcad'])


    def test_misc(self):
        cs = CountingSortExp([5, -3, 12, 0, -3, 101, 15, -23, 51, -44, 2], exp=0)
        res = cs.sort()
        self.assertEqual(res, [-44, -3, -3, -23, 0, 101, 51, 12, 2, 5, 15])

        cs = CountingSortExp([5, -3, 12, 0, -3, 101, 15, -23, 51, -44, 2], exp=0, add_orig_ind=True)
        res = cs.sort()
        self.assertEqual(res, [(-44, 9), (-3, 1), (-3, 4), (-23, 7), (0, 3), (101, 5), (51, 8), (12, 2), (2, 10), (5, 0), (15, 6)])

    def test_second(self):
        cs = CountingSortExp([5, -3, 12, 0, -3, 101, 15, -23, 51, -44, 2], exp=1)
        res = cs.sort()
        self.assertEqual(res, [-44, -23, 5, -3, 0, -3, 101, 2, 12, 15, 51])

    def test_third(self):
        cs = CountingSortExp([5, -3, 12, 0, -3, 101, 15, -23, 51, -44, 2], exp=2)
        res = cs.sort()
        self.assertEqual(res, [5, -3, 12, 0, -3, 15, -23, 51, -44, 2, 101])



class CountingSortExpSpecialListTests(unittest.TestCase):

    def test_string(self):
        with self.assertRaises(Exception) as context:
            cs = CountingSortExpSpecialList('banana')
            res = cs.sort()
        self.assertTrue("Error - Must be list. Type provided is: <class 'str'>" in str(context.exception))

    def test_misc(self):
        with self.assertRaises(Exception) as context:
            cs = CountingSortExpSpecialList([5, -3, 12, 0, -3, 101, 15, -23, 51, -44, 2], exp=0)
            res = cs.sort()
        self.assertTrue("Error: Must contain tuples of integer ranks" in str(context.exception))

        cs = CountingSortExpSpecialList([(11,12,12,11,12), (12,), (11,12), (12,9), (11,)], exp=4, add_orig_ind=True)
        res = cs.sort()
        self.assertEqual(res, [((12,), 1), ((11,12), 2), ((12,9), 3), ((11,), 4), ((11,12,12,11,12), 0)])

        cs = CountingSortExpSpecialList([(11,12,12,11,12), (12,), (11,12), (12,9), (11,)], exp=0, add_orig_ind=True)
        res = cs.sort()
        self.assertEqual(res, [((11,12,12,11,12), 0), ((11,12), 2), ((11,), 4), ((12,), 1), ((12,9), 3)])


class RadixSortTests(unittest.TestCase):

    def test_numeric(self):
        cs = RadixSort([0, 5, -1, 100, 10, -3000, 2])
        res = cs.sort()
        self.assertEqual(res,[-3000, -1, 0, 2, 5, 10, 100])

    def test_string(self):
        cs = RadixSort(['abc', 'a', 'bcae', 'bcad', 'ab'])
        res = cs.sort()
        self.assertEqual(res,['a', 'ab', 'abc', 'bcad', 'bcae'])

    def test_misc(self):
        rs = RadixSort([5, -3, 12, -3000, 0, -3, 101, 15, -23, 51, -44, 2])
        res = rs.sort()
        self.assertEqual(res, [-3000, -44, -23, -3, -3, 0, 2, 5, 12, 15, 51, 101])

        rs = RadixSort([5, -3, 12, -3000, 0, -3, 101, 15, -23, 51, -44, 2], add_orig_ind=True)
        res = rs.sort()
        self.assertEqual(res, [(-3000, 3), (-44, 10), (-23, 8), (-3, 1), (-3, 5), (0, 4), (2, 11), (5, 0), (12, 2), (15, 7), (51, 9), (101, 6)])


        rs = RadixSort(['def', 'adc', 'abc', 'acb', 'acb', 'cab'], all_same_len=True)
        res = rs.sort()
        self.assertEqual(res, ['abc', 'acb', 'acb', 'adc', 'cab', 'def'])

        rs = RadixSort(['def', 'adc', 'abc', 'acb', 'acb', 'cab'], add_orig_ind=True, all_same_len=True)
        res = rs.sort()
        self.assertEqual(res,[('abc', 2), ('acb', 3), ('acb', 4), ('adc', 1), ('cab', 5), ('def', 0)] )

        rs = RadixSort(['ac', 'a', 'ab'])
        res = rs.sort()
        self.assertEqual(res, ['a', 'ab', 'ac'])

        #
        rs = RadixSort(['abc', 'ab', 'a', 'abcd', 'abcb', 'acb'], add_orig_ind=True, all_same_len=False)
        res = rs.sort()
        self.assertEqual(res, [('a', 2), ('ab', 1), ('abc', 0), ('abcb', 4), ('abcd', 3), ('acb', 5)])

        rs = RadixSort(['abc', 'ab', 'a', 'abcd', 'abcb', 'acb'], add_orig_ind=True, all_same_len=True)
        # Expect wrong answer, as not all same length
        res = rs.sort()
        self.assertEqual(res, [('a', 2), ('ab', 1), ('abc', 0), ('abcd', 3),  ('abcb', 4), ('acb', 5)])


class RadixSortSpecialListTests(unittest.TestCase):

    def test_numeric(self):
        with self.assertRaises(Exception) as context:
            rs = RadixSortSpecialList([0, 5, -1, 100, 10, -3000, 2])
            res = rs.sort()
        self.assertTrue("Error: Must be list of rank tuples" in str(context.exception))

    def test_misc(self):
        rs = RadixSortSpecialList([(11,12,12,11,12),(12,), (11,12), (12,9), (11,)], add_orig_ind=True)
        res = rs.sort()
        self.assertEqual(res, [((11,), 4), ((11,12), 2), ((11,12,12,11,12), 0), ((12,), 1), ((12,9), 3)])

        rs = RadixSortSpecialList(
            [(11,12,12,11,12), (14,7), (12,0), (12,9,3), (12,), (11,12), (12,9), (12,-8), (12,9), (12,0), (11,)],
            add_orig_ind=True)
        res = rs.sort()
        self.assertEqual(res, [((11,), 10), ((11,12), 5), ((11,12,12,11,12), 0), ((12,-8), 7), ((12,0), 2), ((12,), 4), ((12,0), 9), ((12,9), 6), ((12,9), 8), ((12,9,3), 3), ((14,7), 1)])

        rs = RadixSortSpecialList([(11, 12, 12, 11, 12), (12,), (11, 12), (12, 9), (11,)], add_orig_ind=False)
        res = rs.sort()
        self.assertEqual(res, [(11,), (11, 12), (11, 12, 12, 11, 12), (12,), (12, 9)])


class SuffixArrayTests(unittest.TestCase):
    def test_misc(self):
        s = 'C'
        sa = SuffixArray(s)
        # was showing key error
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        s = 'CCGTAAACACATGCAGGCAC'
        sa = SuffixArray(s)
        # was showing key error
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

    def test_multi(self):

        for i in range(1000):
            s_len = random.randint(1,100)
            s = gen_rand_string(s_len)
            #print(s)
            sa = SuffixArray(s)
            res = sa.do_it()
            res2 = sa.do_it_naive()
            self.assertEqual(res, res2)

    def test_misc_other(self):

        sa = SuffixArray('AGAGCAAAATTAGCGGAGC')
        # was showing  2,16 when should be 16, 2
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('AGTGAGAGAGCAAAAAAATTAGCGGAGC')
        # 3 nums were wrong - was showing  (8, 20, 25) now  should be (25, 8, 20). Was due to no hash sign in u
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('AGCGACTACTAAGGTAA')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('bananamericanpanamabananas')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('panamabananas')
        # should be [5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12]
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('mississippi')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('yabbadabbado')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArray('banana')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        very_small_sample_s = 'AAAAGTCG'

        sa = SuffixArray(very_small_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        small_sample_s = 'AAAAGTCGCTGCAGCGT'

        sa = SuffixArray(small_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        # correct!  : 0, 1, 2, 12, 3, 11, 6, 14, 8, 10, 13, 7, 15, 4, 16, 5, 9
        self.assertEqual(res, res2)

        slightly_bigger_sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAA'
        # correct!: [27, 26, 0, 1, 2, 19, 21, 23, 12, 3, 25, 11, 17, 6, 14, 8, 18, 20, 22, 24, 10, 13, 7, 15, 4, 16, 5, 9]
        sa = SuffixArray(slightly_bigger_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)


        sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAAAAAAATTAGGCGATGCGGAGC'
        sa = SuffixArray(sample_s)
        # wrong: dc3 is sorting AGCAAAAA before AGC
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        bw_sa_10000_lim = file_to_list('../test_data/large/ecoli_500000.sa')
        bw_sa_dc3 = file_to_list('../test_data/large/ecoli_500000.sa1')
        self.assertEqual(bw_sa_10000_lim[1:], bw_sa_dc3)


class SuffixArrayRankListTests(unittest.TestCase):
    def test_misc(self):
        s = 'C'
        sa = SuffixArrayRankList(s)
        # was showing key error
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        s = 'CCGTAAACACATGCAGGCAC'
        sa = SuffixArrayRankList(s)
        # was showing key error
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

    def test_multi(self):

        for i in range(1000):
            s_len = random.randint(1,100)
            s = gen_rand_string(s_len)
            #print(s)
            sa = SuffixArrayRankList(s)
            res = sa.do_it()
            res2 = sa.do_it_naive()
            self.assertEqual(res, res2)

    def test_misc_other(self):

        sa = SuffixArrayRankList('AGAGCAAAATTAGCGGAGC')
        # was showing  2,16 when should be 16, 2
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('AGTGAGAGAGCAAAAAAATTAGCGGAGC')
        # 3 nums were wrong - was showing  (8, 20, 25) now  should be (25, 8, 20). Was due to no hash sign in u
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('AGCGACTACTAAGGTAA')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('bananamericanpanamabananas')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('panamabananas')
        # should be [5, 3, 1, 7, 9, 11, 6, 4, 2, 8, 10, 0, 12]
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('mississippi')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('yabbadabbado')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        sa = SuffixArrayRankList('banana')
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        very_small_sample_s = 'AAAAGTCG'

        sa = SuffixArrayRankList(very_small_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        small_sample_s = 'AAAAGTCGCTGCAGCGT'

        sa = SuffixArrayRankList(small_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        # correct!  : 0, 1, 2, 12, 3, 11, 6, 14, 8, 10, 13, 7, 15, 4, 16, 5, 9
        self.assertEqual(res, res2)

        slightly_bigger_sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAA'
        # correct!: [27, 26, 0, 1, 2, 19, 21, 23, 12, 3, 25, 11, 17, 6, 14, 8, 18, 20, 22, 24, 10, 13, 7, 15, 4, 16, 5, 9]
        sa = SuffixArrayRankList(slightly_bigger_sample_s)
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)


        sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAAAAAAATTAGGCGATGCGGAGC'
        sa = SuffixArrayRankList(sample_s)
        # wrong: dc3 is sorting AGCAAAAA before AGC
        res = sa.do_it()
        res2 = sa.do_it_naive()
        self.assertEqual(res, res2)

        bw_sa_10000_lim = file_to_list('../test_data/large/ecoli_500000.sa')
        bw_sa_dc3 = file_to_list('../test_data/large/ecoli_500000.sa1')
        self.assertEqual(bw_sa_10000_lim[1:], bw_sa_dc3)


class PopDistTests(unittest.TestCase):

    def test_hw_ideal(self):
        F = 0.0
        init_a = 0.2
        pop = None
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop,  F = F, verbose=1)
        self.assertEqual([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(10)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, (math.pow(1-init_a,2), 2 * init_a * (1 - init_a), math.pow(init_a,2)))
        x = 1

    def test_hw(self):
        F = 0.0
        init_a = 0.2
        pop = 10000
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop, F=F, verbose=1)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(10)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 0.8097, 0.1903])

        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, (math.pow(1-init_a,2), 2 * init_a * (1 - init_a), math.pow(init_a,2)))

    def test_hw_ideal_with_F(self):
        F = 0.4
        init_a = 0.2
        pop = None
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop,  F = F, verbose=1)
        self.assertEqual([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(10)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, [0.704, 0.192, 0.104])
        x = 1

    def test_hw_with_F(self):
        F = 0.4
        init_a = 0.2
        pop = 10000
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop, F=F, verbose=1)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa],
                                       [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(10)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa],
                                       [1 - init_a, init_a, 0.8055, 0.1945])

        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, [0.704, 0.192, 0.104])

    def test_fitness_aa_ideal(self):
        F = 0.0
        init_a = 0.2
        pop = None
        fitnesses = [1.0, 1.0, 0.9]
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop,  F = F,  genotype_fitnesses=fitnesses, verbose=1)
        self.assertEqual([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(10)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 0.828441, 0.171559])
        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, (math.pow(1-init_a,2), 2 * init_a * (1 - init_a), math.pow(init_a,2)))
        x = 1

    def test_very_fit_het(self):
        F = 0.0
        init_a = 0.2
        pop = None
        fitnesses = [0.9, 1.0, 0.9]
        np.random.seed(42)
        pdhw = PopDist(init_a, pop=pop,  F = F,  genotype_fitnesses=fitnesses, verbose=1)
        self.assertEqual([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 1 - init_a, init_a])
        pdhw.sim_generations(100)
        np.testing.assert_almost_equal([pdhw.init_fA, pdhw.init_fa, pdhw.current_fA, pdhw.current_fa], [1 - init_a, init_a, 0.5016925, 0.4983075])
        res = PopDistGen.hw_genotype_freqs(init_a, F=F)
        np.testing.assert_almost_equal(res, (math.pow(1-init_a,2), 2 * init_a * (1 - init_a), math.pow(init_a,2)))
        x = 1

    def test_json_dump(self):
        np.random.seed(42)
        pd = PopDist(0.5, genotype_fitnesses=[0.9, 1.0, 1.0], pop=100, F=0, verbose=1)
        print(pd)
        pd.sim_generations(3)
        j_out = pd.to_json()
        self.assertEqual(j_out, '{"init_fa": 0.5, "init_pop": 100, "init_genotype_fitnesses": [0.9, 1.0, 1.0], "init_F": 0, "verbose": 1, "gens": [{"in_fa": 0.5, "pop": 100, "genotype_fitnesses": [0.9, 1.0, 1.0], "F": 0, "gen_num": 0, "prev_gen_num": null, "verbose": 1, "random_mating_genotypes": null, "out_fa": 0.5}, {"in_fa": 0.5, "pop": 100, "genotype_fitnesses": [0.9, 1.0, 1.0], "F": 0, "gen_num": 1, "prev_gen_num": 0, "verbose": 1, "random_mating_genotypes": null, "random_mating_genotype_counts": [24, 44, 32], "random_mating_genotype_freqs": [0.24, 0.44, 0.32], "survived_genotype_freqs": [0.22131147540983603, 0.4508196721311475, 0.32786885245901637], "out_fa": 0.5532786885245902}, {"in_fa": 0.5532786885245902, "pop": 100, "genotype_fitnesses": [0.9, 1.0, 1.0], "F": 0, "gen_num": 2, "prev_gen_num": 1, "verbose": 1, "random_mating_genotypes": null, "random_mating_genotype_counts": [22, 47, 31], "random_mating_genotype_freqs": [0.22, 0.47, 0.31], "survived_genotype_freqs": [0.20245398773006137, 0.4805725971370143, 0.3169734151329243], "out_fa": 0.5572597137014315}, {"in_fa": 0.5572597137014315, "pop": 100, "genotype_fitnesses": [0.9, 1.0, 1.0], "F": 0, "gen_num": 3, "prev_gen_num": 2, "verbose": 1, "random_mating_genotypes": null, "random_mating_genotype_counts": [16, 50, 34], "random_mating_genotype_freqs": [0.16, 0.5, 0.34], "survived_genotype_freqs": [0.14634146341463417, 0.508130081300813, 0.3455284552845529], "out_fa": 0.5995934959349594}]}')
        pd_inflated = PopDist.pop_dist_from_json(j_out)
        print('inf: ', pd_inflated)
        np.testing.assert_almost_equal([pd.current_fA, pd.current_fa], [pd_inflated.current_fA, pd_inflated.current_fa])

class CrossTests(unittest.TestCase):

    def testsimple(self):
        random.seed(42)
        c1 = ChromosomeTemplate.from_symbol_list(['A-R-1000', 'B-D-2000', 'C-R-100000000', 'D-R-200000000'], 'Chrom 1')
        gt = GenomeTemplate(ploidy=2, chromosome_templates=[c1], name='Onechrom')
        parent_1 = Organism.organism_with_het_genotype(gt)
        parent_2 = Organism.organism_with_het_genotype(gt)
        last_child_added = parent_1.mate(parent_2, times=4)
        self.assertEqual(str(parent_1.genome), 'AaBbCcDd')
        self.assertEqual(str(parent_2.genome), 'AaBbCcDd')
        self.assertEqual(parent_1.num_children, 4)
        self.assertEqual(parent_2.num_children, 4)
        self.assertEqual(last_child_added.num_siblings, 3)
        self.assertEqual(str(parent_1.children[0].genome), 'AABBCcDd')
        self.assertEqual(str(parent_1.children[1].genome), 'AABBCCDd')
        self.assertEqual(str(parent_1.children[2].genome), 'aabbCcDD')
        self.assertEqual(str(parent_1.children[3].genome), 'aabbCcDd')
        x = 1


    def test_two_chroms(self):
        random.seed(42)
        c1 = ChromosomeTemplate.from_symbol_list(['A-R-1000', 'B-D-2000'], 'Chrom 1 Linked Genes')
        c2 = ChromosomeTemplate.from_symbol_list(['C'], 'Chrom 2')
        gt = GenomeTemplate(ploidy=2, chromosome_templates=[c1, c2], name='Twochroms')
        parent_1 = Organism.organism_with_het_genotype(gt)
        parent_2 = Organism.organism_with_het_genotype(gt)
        last_child_added = parent_1.mate(parent_2, times=8)
        self.assertEqual(str(parent_1.genome), 'AaBbCc')
        self.assertEqual(str(parent_2.genome), 'AaBbCc')
        self.assertEqual(parent_1.num_children, 8)
        self.assertEqual(parent_2.num_children, 8)
        self.assertEqual(last_child_added.num_siblings, 7)
        self.assertEqual(str(parent_1.children[0].genome), 'AABBCC')
        self.assertEqual(str(parent_1.children[1].genome), 'AABBCC')
        self.assertEqual(str(parent_1.children[2].genome), 'AaBbCc')
        self.assertEqual(str(parent_1.children[3].genome), 'AABBcc')
        x = 1


if __name__ == '__main__':
    unittest.main()

