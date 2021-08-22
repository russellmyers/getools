import numpy as np
import json

class Sequence:

    def __init__(self,seq, name = ''):
        self.seq = seq
        self.name = name
        self.alphabet = None

    def __str__(self):
        return(self.seq)

    def length(self):
        return len(self.seq)

    def is_valid(self):
        return True


    def num_diffs(self, other_seq):
        diffs = 0
        for i, ch in enumerate(self.seq):
            if other_seq.seq[i] != ch:
                diffs +=1
        return diffs



class  DNASequence(Sequence):

    DNA_ALPHABET  = ['A','C','G','T']

    SIMPLE_SCORING_MATRIX = [[1,-2,-1,-2],
                              [-2,1,-2,-1],
                              [-1,-2,1,-2],
                              [-2,-1,-2,1]]


    def __init__(self,seq):
        super().__init__(seq)
        self.alphabet = DNASequence.DNA_ALPHABET

    def is_valid(self):
        for base in self.seq:
            if base in self.alphabet:
                pass
            else:
                return False
        return True


class RNASequence(Sequence):
    RNA_ALPHABET = ['A','C','G','U']

    def __init__(self, seq):
        super().__init__(seq)
        self.alphabet = RNASequence.RNA_ALPHABET

    def is_valid(self):
        for base in self.seq:
            if not base in self.alphabet:
                return False
        return True

class ScoringMatrix:

    def __init__(self, mat, alphabet):
        self.alphabet = alphabet
        self.mat = mat
        self.inds = {ch:pos for pos, ch in enumerate(self.alphabet)}

    def score(self, c1, c2):
       return self.mat[self.inds[c1]][self.inds[c2]]


class Aligner:

    INDEL_CHAR = '-'

    def __init__(self,scoring_matrix=None,mismatch_pen=-1,indel_pen=-1, match_score = 1, verbose=0):
        self.s = None
        self.t = None
        self.matrix = None
        self.scoring_matrix = scoring_matrix
        self.mismatch_pen = mismatch_pen
        self.indel_pen = indel_pen
        self.match_score = match_score
        self.verbose = verbose

    def _init_matrix(self):
        matrix = np.zeros((self.s.length() + 1, self.t.length() + 1), dtype=int)
        if self.verbose > 0:
            print(matrix.shape)
        return matrix

    def _score_at_pos(self, i, j):
        return self.matrix[i][j]


    def calc_match_score(self,c1,c2):
        if self.scoring_matrix is None:
           if c1 == c2:
               return self.match_score
           else:
               return self.mismatch_pen
        else:
           #i = self.s.alphabet.index(c1)
           #j = self.t.alphabet.index(c2)
           #return self.scoring_matrix[i][j]
           return self.scoring_matrix.score(c1, c2)


    def _is_match(self, i, j):
        if (i == 0) or (j == 0):
            return False

        return str(self.s)[i - 1] == str(self.t)[j - 1]

    def _build_scores(self):
        return

    def backtrack(self):
        pass

    def align(self, s, t):
        self.s = s
        self.t = t
        self.matrix = self._init_matrix()
        self._build_scores()
        return self.matrix[-1][-1]

    def add_variant(self,variants,pos,ref,alt):
        variants.append([pos+1,ref,alt])

    def get_next_var(self,ref,t, prev_ref_base, prev_t_base):
        if ref[0] == Aligner.INDEL_CHAR:
           ref_base = ref[0]
           indel_length = 0
           while (indel_length < len(ref)-1) and (ref_base == Aligner.INDEL_CHAR):
               indel_length +=1
               ref_base = ref[indel_length]
           if indel_length >=  len(ref) - 1:
               indel_length = len(ref)
           if prev_ref_base is None:
              return(indel_length,ref[indel_length], t[:indel_length+1])
           else:
              return (indel_length, prev_ref_base, prev_ref_base + t[:indel_length])
        elif t[0] == Aligner.INDEL_CHAR:
            t_base = t[0]
            indel_length = 0
            while (indel_length < len(ref)-1) and (t_base == Aligner.INDEL_CHAR):
                indel_length += 1
                t_base = t[indel_length]
            if indel_length >= len(ref) - 1:
                indel_length = len(ref)
            if prev_t_base is None:
                return (indel_length, ref[:indel_length+1], t[indel_length])
            else:
                return (indel_length, prev_ref_base + ref[:indel_length], prev_ref_base)


        elif ref[0] != t[0]:
               return 1,ref[0],t[0]
        else:
            return 1,None,None

    def al_to_vcf(self,al_ref,al_t):
        variants = []
        curr_pos = 0
        while curr_pos < len(al_ref):
            offset,ref,alt = self.get_next_var(al_ref[curr_pos:],al_t[curr_pos:], al_ref[curr_pos - 1] if curr_pos > 0 else None,al_t[curr_pos - 1] if curr_pos > 0 else None)
            if ref is None:
                pass
            else:
                self.add_variant(variants,curr_pos,ref,alt)
            curr_pos += offset

        return variants


    # def alignment_to_vcf(self,al_ref,al_t):
    #     variants = []
    #     done = False
    #     pos = 0
    #     for pos in range(0,len(al_ref)):
    #           ref_base = al_ref[pos]
    #           t_base = al_t[pos]
    #
    #           if ref_base == t_base:
    #              continue
    #
    #           if ref_base == Aligner.INDEL_CHAR: # Insertion in t
    #              if pos == 0:
    #                 ref = al_ref[pos+1]
    #                 alt = al_t[pos]
    #              else:
    #                 ref = al_ref[pos-1]
    #                 alt = al_t[pos-1]
    #
    #              offset = 1
    #              while ref_base == Aligner.INDEL_CHAR:
    #                  alt += al_t[pos + offset]
    #                  offset+=1
    #                  if (pos + offset) >= len(al_ref):
    #                      ref_base = '*'
    #                  else:
    #                     ref_base = al_ref[pos + offset]
    #           elif t_base == Aligner.INDEL_CHAR: # Deletion in t
    #              if pos == 0:
    #                 ref = al_ref[pos]
    #                 alt = al_t[pos+1]
    #              else:
    #                 ref = al_ref[pos]
    #                 alt = al_t[pos]
    #              offset = 1
    #              while t_base == Aligner.INDEL_CHAR:
    #                  ref += al_ref[pos + offset]
    #                  offset += 1
    #                  if (pos + offset)  >= len(al_ref):
    #                      t_base = '*'
    #                  else:
    #                     t_base = al_t[pos + offset]
    #           elif ref_base != t_base:
    #              ref = al_ref[pos]
    #              alt = al_t[pos]
    #
    #
    #           self.add_variant(variants,pos,ref,alt)
    #
    #     return variants



    def __str__(self):
        return (f'Aligner - mismatch pen: {self.mismatch_pen} indel pen: {self.indel_pen}')


class GlobalAligner(Aligner):



    def _build_scores(self):

        if self.verbose > 0:
           print('Building scores:')
        for i in range(self.matrix.shape[0]):
            if i % 100 == 0:
                if self.verbose > 0:
                    print('Processing row: ',i, '/',self.matrix.shape[0], i / self.matrix.shape[0] * 100, '%')
            for j in range(self.matrix.shape[1]):
                if i == 0 and j == 0:
                    continue
                elif i == 0:
                    self.matrix[i][j] = self._score_at_pos(i,j-1) + self.indel_pen
                elif j == 0:
                    self.matrix[i][j] = self._score_at_pos(i-1, j) + self.indel_pen
                else:
                    self.matrix[i][j] = max(
                        self._score_at_pos(i-1,j-1) + (self.calc_match_score(str(self.s)[i-1],str(self.t)[j-1])),
                        self._score_at_pos(i - 1, j) + self.indel_pen,
                        self._score_at_pos(i,j-1) + self.indel_pen
                    )

    def backtrack(self):
        al_s = []
        al_t = []
        i = self.matrix.shape[0] - 1
        j = self.matrix.shape[1] -1
        done = False
        while not done:
            if i == 0 and j == 0:
                done = True
            elif i == 0:
                al_s.insert(0,Aligner.INDEL_CHAR)
                al_t.insert(0,str(self.t)[j-1])
                j -=1
            elif j == 0:
                al_s.insert(0,str(self.s)[i - 1])
                al_t.insert(0,  Aligner.INDEL_CHAR)
                i -= 1
            else:
                if self._score_at_pos(i-1,j-1) + self.calc_match_score(str(self.s)[i-1],str(self.t)[j-1]) == self._score_at_pos(i,j):
                    al_s.insert(0,str(self.s)[i-1])
                    al_t.insert(0,str(self.t)[j-1])
                    i -=1
                    j -=1
                elif  self._score_at_pos(i-1,j) + self.indel_pen == self._score_at_pos(i,j):
                    al_s.insert(0, str(self.s)[i-1])
                    al_t.insert(0, Aligner.INDEL_CHAR)
                    i -=1
                elif self._score_at_pos(i, j-1) + self.indel_pen == self._score_at_pos(i, j):
                    al_s.insert(0, Aligner.INDEL_CHAR)
                    al_t.insert(0, str(self.t)[j-1])
                    j -=1
                else:
                    raise(Exception('Unable to determine backtrack at position: ' + str(i) + ' ' + str(j)))


        return [''.join(al_s),''.join(al_t)]



class QuickAligner:
    '''
    Class to check if sequences are basically the same barring SNPs (hence can be compared base by base for SNPS) but start at different positions

    '''

    def __init__(self, verbose=0):
        self.verbose = verbose

    def num_indels_prefix(self,s):
        num_indels_at_start = 0
        for ch in s:
            if ch == Aligner.INDEL_CHAR:
                num_indels_at_start +=1

        return num_indels_at_start

    def start_matching_part(self,al_s,al_t):

        full_match_start = 0

        for i in range(len(al_s)-200,-1,-1):
            if (al_s[i] == Aligner.INDEL_CHAR) or (al_t[i] == Aligner.INDEL_CHAR):
                full_match_start = i + 1
                break

        s_offset = full_match_start - self.num_indels_prefix(al_s[:full_match_start])
        t_offset = full_match_start - self.num_indels_prefix(al_t[:full_match_start])

        return s_offset, t_offset


    def quick_align(self,s,t,s_offset,t_offset):
        s1 = str(s)[s_offset:]
        t1 = str(t)[t_offset:]
        n_positions = []
        diffs = []
        if self.verbose > 0:
            print('comparing from: ', s1[:50], t1[:50])
        patience = None
        equal_run = 0
        for i in range(max(len(s1), len(t1))):
            actual = i + s_offset
            if patience is None:
                pass
            else:
                patience += 1

            if i >= len(t1):
                if self.verbose > 0:
                    print('diff extra seq1: ', actual, s1[i:])
                diffs.append({'pos':actual,'type':'extra seq1','ref':s1[i:],'other':''})
                break
            if i >= len(s1):
                if self.verbose > 0:
                    print('diff extra seq2: ', actual, t1[i:])
                diffs.append({'pos': actual, 'type': 'extra seq2', 'ref': '', 'other': t1[i:]})
                break
            if t1[i] == 'N':
                n_positions.append(actual)
                equal_run += 1  # don't reset for this
            elif s1[i] == t1[i]:
                equal_run += 1
            else:
                if self.verbose > 0:
                    print('diff: ', actual, s1[i], t1[i])
                diffs.append({'pos':actual,'type':'SNP','ref':s1[i],'other':t1[i]})
                if patience is None:
                    patience = 1  # start counting
                    equal_run = 0
                else:
                    if patience > 30:
                        if equal_run < 10:
                            print('Failed - equal run was: ', equal_run)
                            return None
                            break
                        else:
                            patience = 1  # start counting again
                            equal_run = 0
                    else:
                        equal_run = 0
        return diffs

    def align(self,s,t, prealign_size = 500, equal_thresh=50):

        diffs = None

        al = GlobalAligner(scoring_matrix=s.SIMPLE_SCORING_MATRIX, verbose=0, indel_pen=-3)
        s1 = DNASequence(str(s)[:prealign_size])
        t1 = DNASequence(str(t)[:prealign_size])

        al_result = al.align(s1, t1)
        if self.verbose > 0:
            print(al_result)

        al_s, al_t = al.backtrack()
        if self.verbose > 0:
            print(al_s)
            print(al_t)


        #s_offset = self.num_indels_prefix(al_t)
        #t_offset = self.num_indels_prefix(al_s)
        s_offset,t_offset = self.start_matching_part(al_s, al_t)

        if self.verbose > 0:
            print(s_offset,t_offset)
        num_bases_equal = 0
        for i in range(equal_thresh):
            if str(s)[s_offset + i] == str(t)[t_offset + i]:
                num_bases_equal +=1
        if num_bases_equal >= equal_thresh - 10:
            print('Candidate for quick alignment')
            diffs = self.quick_align(s, t, s_offset, t_offset)
        else:
            print('Not candidate for quick alignment',str(s)[s_offset:s_offset + equal_thresh],str(t)[t_offset:t_offset + equal_thresh] )

        return diffs

if __name__  == '__main__':

    seq1 = Sequence('panamabananas')
    seq2 = Sequence('anamaqqbananas')
    print(seq1, seq1.is_valid())
    al = GlobalAligner()
    score = al.align(seq1, seq2)
    res = al.backtrack()


    seq1 = DNASequence('AGCTTAGCT')
    seq2 = DNASequence('ACCGAGCT')
    print(seq1, seq1.is_valid())
    al = GlobalAligner(scoring_matrix=ScoringMatrix(DNASequence.SIMPLE_SCORING_MATRIX, DNASequence.DNA_ALPHABET))
    score = al.align(seq1, seq2)
    res = al.backtrack()

    seq = DNASequence('CGACTGTGTAGCTGATGHCTGTAGTCGTAGCTGATCGTACG')
    print(seq, seq.is_valid())

    s = DNASequence('TACGTTTTTTTTCCCCCAAAAA')
    t = DNASequence('ACGTTTTAGGGATTTCCAAAAA')

    with open("../test_data/covid19/1_ref_sequence_NC_045512.2.fasta") as file:  # Use file to refer to the file object
        data = file.read()
        data_list = data.split('\n')
        seq_1 = ''.join(data_list[1:])
        print('Length ref sequence: ',len(seq_1), ' First 50: ',seq_1[:50])

    import os

    use_all = False

    if use_all:
        sequence_files = []
        for file in os.listdir("../test_data/covid19"):
            if file.endswith(".fasta"):
               file_num = file.split('_')[0]
               if file_num == '1':
                   pass # ref
               else:
                   sequence_files.append(file)
    else:
        # sequence_files = ['12_MN996532.1.fasta']
        # name_abbrev = 'vars_ratg13.json'
        sequence_files = ['14_MK334047.1.fasta']
        name_abbrev = 'vars_human_coronavirus_nl63.json'

    print(sequence_files)


    all_diffs = []
    for file_name in sequence_files:
        with open("../test_data/covid19/" + file_name) as file:  # Use file to refer to the file object
            print('Processing: ',file_name)
            data = file.read()
            data_list = data.split('\n')
            seq_2 = ''.join(data_list[1:])
            #print('Length comparison sequence: ', len(seq_2), ' First 50: ', seq_2[:50])

            if use_all:
                qa = QuickAligner(verbose=0)
                diffs = qa.align(DNASequence(seq_1), DNASequence(seq_2))
                if diffs is None:
                    print('diffs not calced')
                else:
                    print('diffs: ', diffs)
                    all_diffs.append({'file': file_name, 'diffs': diffs})

                with open("../test_data/covid19/" + "out_diffs.txt",
                          'w') as file:  # Use file to refer to the file object
                    for diff_rec in all_diffs:
                        file.write(diff_rec['file'] + '\t' + json.dumps(diff_rec['diffs']) + '\n')

            else:
                qa = GlobalAligner(verbose=1)
                seq_1 = seq_1[:5000]
                seq_2 = seq_2[:5000]
                print('gl align: ',qa.align(DNASequence(seq_1), DNASequence(seq_2)))
                al_s, al_t = qa.backtrack()
                print(al_s[:100])
                print(al_t[:100])
                variants = qa.al_to_vcf(al_s, al_t)
                print('len vars: ', len(variants))
                with open("../test_data/covid19/" + name_abbrev, "w") as write_file:
                    json.dump(variants, write_file)




    exit()



    #file_name = '6_MN938384.1.fasta' # '11_MT118835.1.fasta' #'6_MN938384.1.fasta' #'5_LC528233.1.fasta' #'2_MN985325.1.fasta'  # '3_MT246490.1.fasta' # '4_MT233522.1.fasta'

    with open("../test_data/covid19/" + file_name) as file:  # Use file to refer to the file object
        data = file.read()
        data_list = data.split('\n')
        seq_2 = ''.join(data_list[1:])
        print('Length comparison sequence: ',len(seq_2), ' First 50: ',seq_2[:50])


    #s = DNASequence('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTG')
    #t = DNASequence('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTG')

    s = DNASequence(seq_1)
    t = DNASequence(seq_2)

    qa = QuickAligner(verbose = 1)
    diffs = qa.align(s,t)

    ref_start_offset = 0
    t_start_offset = 0
    s1 = seq_1[ref_start_offset:]
    s2 = seq_2[t_start_offset:]

    n_positions = []
    print('comparing from: ',s1[:50], s2[:50])
    patience = None
    equal_run = 0
    for i in range(max(len(s1),len(s2))):
        actual = i + ref_start_offset
        if patience is None:
            pass
        else:
            patience +=1

        if i >= len(s2):
           print('diff extra seq1: ',actual,s1[i:])
           break
        if i >= len(s1):
            print('diff extra seq2: ', actual, s2[i:])
            break
        if s2[i] == 'N':
           n_positions.append(actual)
           equal_run +=1 #don't reset for this
        elif s1[i] == s2[i]:
            equal_run +=1
        else:
            print('diff: ',actual,s1[i],s2[i])
            if patience is None:
                patience = 1 # start counting
                equal_run = 0
            else:
                if patience > 30:
                    if equal_run < 10:
                       print('equal run was: ',equal_run)
                       break
                    else:
                       patience = 1 # start counting again
                       equal_run = 0
                else:
                     equal_run = 0


    print('i: ',i)
    s = DNASequence(seq_1[:500])
    t = DNASequence(seq_2[:500])


    al = GlobalAligner(scoring_matrix = s.SIMPLE_SCORING_MATRIX, verbose = 1, indel_pen=-3)

    print(al)
    print(al.align(s,t))

    al_s,al_t = al.backtrack()
    print(al_s)
    print(al_t)

    variants = al.al_to_vcf(al_s,al_t)
    print('vars: ',variants)
    x=1