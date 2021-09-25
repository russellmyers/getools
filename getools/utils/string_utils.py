""" Various string utilities - fast searching algorithms such as  etc"""

import random
import time
from getools.utils.file_utils import read_fasta, string_to_file, read_txt_file, list_to_file, file_to_list
from collections import defaultdict
import math


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print(f'{method.__name__}  {(te - ts) :.4f} secs')
        return result
    return timed


@timeit
def pattern_match_naive(s, t):
    matches = []
    iters = 0
    for i in range(len(s) - len(t) + 1):
        match = True
        for j, t_ch in enumerate(t):
            iters += 1
            if t_ch != s[i+j]:
                match = False
                break
        if match:
            matches.append(i)
    print(f'Iters actually performed: {iters}')
    return matches


@timeit
def pattern_match_zbox(s, t=None, also_return_zbox=False, print_iters=False):

    _l = -1
    r = -1
    _Z = []

    matches = []
    if t is None:
        t = ''
        c = s
    else:
        c = t + '$' + s
    len_c = len(c)
    len_t = len(t)

    _Z.append(-1)
    for i in range(1, len_c):
        if i > r:
            c_i = i
            t_i = 0
            num_matches = 0
            while c_i < len_c and c[c_i] == c[t_i]:
                num_matches += 1
                c_i += 1
                t_i += 1
            if num_matches > 0:
                _l = i
                r = i + num_matches - 1
            _Z.append(num_matches)
        else:
            beta = r - i + 1
            rel = i - _l
            if _Z[rel] < beta:
                _Z.append(_Z[rel])
            elif _Z[rel] > beta:
                _Z.append(beta)
                _l = i
            else:
                c_i = r + 1
                t_i = beta
                num_matches = 0
                while c_i < len_c and c[c_i] == c[t_i]:
                    num_matches += 1
                    c_i += 1
                    t_i += 1

                _l = i
                if num_matches > 0:
                    r = r + num_matches
                _Z.append(r - _l + 1)

        if _Z[i] == len(t):
            matches.append(i - len_t - 1)

    if also_return_zbox:
        return matches, _Z
    else:
        return matches


class BWT:
    def __init__(self, s, bw_s=None, sa=None, use_linear_sa_construction=True, max_sort_len_for_naive=None, verbose=0):
        self.verbose = verbose
        if s is None:
            self.s = None  # Reconstruct from bw_s if supplied?
        else:
            self._s = s + '$'
            if self.verbose > 0:
                print(f'BWT for {s[:5]}... Total len: {len(s)}')

        self.max_sort_len_for_naive = max_sort_len_for_naive  # Used in naive sa construction to limit sort lengths. If None - no limit

        if bw_s is None:
            if use_linear_sa_construction:
                self.make_bwt()
            else:
                self.make_bwt_naive_lambda_sort(max_sort_len=self.max_sort_len_for_naive)
        else:
            if self.verbose > 0:
                print('  Already provided bwt and sa')
            self.bw = bw_s
            self.sa = sa
            self._left, self._right = self._calc_first_lasts()

        # # TODO: Replace this with partial suffix array stuff
        # if self.also_return_match_inds:
        #     self.s_prime, self.s_prime_inds = self.reconstruct_naive(also_return_left_inds=True)

    @staticmethod
    def from_bw(bw_s, sa=None, verbose=0):
        return BWT(s=None, bw_s=bw_s, sa=sa, verbose=verbose)

    # @timeit
    # def make_bwt_naive(self):
    #     rotations = []
    #     for i in range(len(self._s), 0, -1):
    #       rotations.append(self._s[i:] + self._s[:i])
    #     rotations.sort()
    #     self._rotations = rotations
    #     self.bw = ''.join([x[-1] for x in rotations])
    #     self._left, self._right = self._calc_first_lasts()

    @timeit
    def make_bwt(self):
        if self.verbose > 0:
            print('  Creating suffix array linear')
        self.sa = SuffixArrayRankList(self._s).do_it()
        if self.verbose > 0:
            print('  Constructing bwt')
        # self.bw = ''.join([(self._s[x:] + self._s[:x])[-1] for x in self.sa])
        self.bw = ''.join(['$' if x == 0 else self._s[x-1] for x in self.sa])
        # for i in range(len(self._s), 0, -1):
        #     rotations.append(self._s[i:] + self._s[:i])
        # rotations.sort()
        # self._rotations = rotations
        # self.bw = ''.join([x[-1] for x in rotations])
        self._left, self._right = self._calc_first_lasts()

    @timeit
    def make_bwt_naive_lambda_sort(self, max_sort_len=None):
        if self.verbose > 0:
            print('  Creating suffix array naive')
        if max_sort_len is None:
            max_sort_len = len(self._s)
        rotations = []
        rotations = [i for i in range(len(self._s))]
        rotations.sort(key=lambda x: (self._s[x:] + self._s[:x])[:max_sort_len])
        self.sa = rotations
        if self.verbose > 0:
            print('  Constructing bwt')
        self.bw = ''.join([(self._s[x:] + self._s[:x])[-1] for x in rotations])
        # for i in range(len(self._s), 0, -1):
        #     rotations.append(self._s[i:] + self._s[:i])
        # rotations.sort()
        # self._rotations = rotations
        # self.bw = ''.join([x[-1] for x in rotations])
        self._left, self._right = self._calc_first_lasts()

    # noinspection PyMethodMayBeStatic
    def save(self, file_name, also_save_suffix_array=False):
        string_to_file(f'{file_name}', bw.bw)
        if also_save_suffix_array:
            ext = file_name[-3:]
            if ext == '.bw':
                list_to_file(f'{file_name[:-3]}.sa', bw.sa)
            else:
                list_to_file(f'{file_name}.sa', bw.sa)

    def _calc_first_lasts(self):
        if self.verbose > 0:
            print('  Calcing first/lasts')
        counts = {}
        count_array = []
        for i, ch in enumerate(self.bw):
            if ch in counts:
                counts[ch] += 1
            else:
                counts[ch] = 1
            count_entry = {key: val for key, val in counts.items()}
            count_array.append(count_entry)

        left = {}
        tot = 0
        sorted_keys = sorted(counts.keys())
        for key in sorted_keys:
            left[key] = {'top': tot, 'bottom': tot + counts[key] - 1}
            tot += counts[key]

        right = count_array
        return left, right

    def find_left_pos(self, ch, num):
        return self._left[ch]['top'] + num

    @timeit
    def reconstruct_naive(self, also_return_suffix_array=False):

        orig = ''
        if also_return_suffix_array:
            suffix_array = [-1 for _ in range(len(self.bw))]
        ind = self._left['$']['top']
        done = False
        while not done:
            r_ch = self.bw[ind]
            if r_ch == '$':
                done = True
            else:
                orig = r_ch + orig
                ind = self.find_left_pos(r_ch, self._right[ind][r_ch]-1)
                if also_return_suffix_array:
                    suffix_array[ind] = len(self.bw)-1 - len(orig)
        if also_return_suffix_array:
            return orig, suffix_array
        else:
            return orig

    @timeit
    def pattern_match(self, t, force_reconstruct=False):
        # if suffix array not available:
        #    - if force_reconstruct is True:  will reconstruct the original string (and suffix array)
        #    - if force_reconstruct is False:will just return bwt index for pattern matches, not actual string positions

        t_ind = len(t) - 1
        l_ch = t[t_ind]
        if l_ch in self._left:
            l_range_top = self._left[l_ch]['top']
            l_range_bottom = self._left[l_ch]['bottom']
        else:
            return []

        done = False
        while not done:
            r_ch = t[t_ind-1]
            if r_ch not in self._left:
                # char not in alphabet
                return []

            if r_ch not in self._right[l_range_top]:
                first_right = 1
            else:
                first_right = self._right[l_range_top][r_ch] if self.bw[l_range_top] == r_ch\
                    else self._right[l_range_top][r_ch] + 1
            if r_ch not in self._right[l_range_bottom]:
                last_right = 0
            else:
                last_right = self._right[l_range_bottom][r_ch]
            l_range_top = self.find_left_pos(r_ch, first_right-1)
            l_range_bottom = self.find_left_pos(r_ch, last_right-1)

            if l_range_top > l_range_bottom:
                return []

            t_ind -= 1
            if t_ind < 1:
                done = True

        candidates = [i for i in range(l_range_top, l_range_bottom+1)]

        if force_reconstruct and self.sa is None:
            print('No suffix array available. Will reconstruct before pattern match')
            _, self.sa = self.reconstruct_naive(also_return_suffix_array=True)

        if self.sa is not None:
            # TODO  use partial suffix array instead
            candidate_inds = [self.sa[c] for c in candidates]
            return candidate_inds
        else:
            # Suffix array not available, just return bwt matches, not char positions
            print('Warning: No suffix array available - returning bwt indexes only.'
                  ' If char positions rqd, add "force_reconstruct=True" param')
            return [f'bwtind:{c}' for c in candidates]

    def compress_bw(self):
        out = ''
        last_ch = ''
        count = 0
        for ch in self.bw:
            if ch == last_ch:
                count += 1
            else:
                if count > 1:
                    out += f'{count}{last_ch}'
                elif count == 1:
                    out += last_ch
                count = 1
            last_ch = ch
        out += f'{count}{last_ch}'
        return out


class SuffixArray:

    def assign_lex_names(self, triplets):

        char_codes = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz{}'

        curr_name = 0
        prev_val = '***'
        triplet_names = {}
        triplet_inds = {}
        all_unique = True

        for i, tup in enumerate(triplets):
            if tup[0] == prev_val:
                all_unique = False
            else:
                curr_name += 1
            triplet_names[tup] = char_codes[curr_name]
            triplet_inds[i] = tup
            prev_val = tup[0]

        return triplet_names, triplet_inds, all_unique

    def make_u(self, triplets_with_names):

        u = ''
        mod_1_inds = [i for i in range(1, len(self._s), 3)]
        mod_2_inds = [i for i in range(2, len(self._s), 3)]
        s_with_sentinel = self._s+'$$'
        for ind in mod_1_inds:
            u += str(triplets_with_names[(s_with_sentinel[ind:ind+3], ind)])
        u += '#'
        for ind in mod_2_inds:
            u += str(triplets_with_names[(s_with_sentinel[ind:ind+3], ind)])

        return u, mod_1_inds + [-1] + mod_2_inds

    def construct_1_2(self):

        s_with_sentinel = (self._s + '$') if len(self._s) % 3 == 1 else (self._s + '$$')   # + '$$$'
        len_s = len(self._s)
        triplets_1_2 = []
        for i in range(len_s):  # (len_s - 2):
            if i % 3 == 0:
                pass
            else:
                triplet = s_with_sentinel[i:i+3]
                triplets_1_2.append((triplet, i))
        triplets_1_2 = RadixSort(triplets_1_2).sort()

        triplet_1_2_names, triplet_1_2_inds, all_unique = self.assign_lex_names(triplets_1_2)

        return triplet_1_2_names, triplet_1_2_inds, all_unique

    # noinspection PyMethodMayBeStatic
    def small_sort(self, u):
        sa = [i for i in range(len(u))]
        sa.sort(key=lambda x: u[x:] + u[:x])
        return sa

    # noinspection PyMethodMayBeStatic
    def reconstruct_from_sorted_u(self, u_sorted, u_inds, triplet_1_2_names, triplet_1_2_inds):
        reconstructed = []
        for i, item in enumerate(u_sorted):
            reconstructed.append(u_inds[item])

        return reconstructed

    # noinspection PyMethodMayBeStatic
    def reconstruct_from_inds(self, u, u_inds, triplet_1_2_inds):
        reconstructed = []
        for key in sorted(triplet_1_2_inds.keys()):
            reconstructed.append(triplet_1_2_inds[key][1])

        return reconstructed

    def get_ranks_3(self, u, sa_12):
        ranks_1_2 = {el: i+1 for i, el in enumerate(sa_12)}
        # if len(u) % 3 == 0:
        #     ranks_3 = [(u[i] + str(ranks_1_2[i + 1]).zfill(3), i) for i in range(0, len(u)-1, 3)]
        #     # TODO - not sure about next line
        #     ranks_3.append((u[-1] + '000', len(ranks_3)))
        # else:
        ranks_3 = [((u[i] + '000'), i) if i >= len(u)-1
                   else (u[i] + str(ranks_1_2[i+1]).zfill(3), i) for i in range(0, len(u), 3)]
        ranks_3 = RadixSort(ranks_3).sort()
        ranks_3_sa = [x[1] for x in ranks_3]
        return ranks_3_sa

    # noinspection PyMethodMayBeStatic,PyBroadException
    def merge(self, s_with_dollar, sa_3, sa_1_2):
        sa_1_2_ranks = {el: i+1 for i, el in enumerate(sa_1_2)}

        merged = []
        ind_3 = 0
        ind_1_2 = 0
        done = False
        # s_with_dollar = self._s + '$$$'
        while not done:
            if ind_1_2 >= len(sa_1_2):
                merged.append(sa_3[ind_3])
                ind_3 += 1
            elif ind_3 >= len(sa_3):
                merged.append(sa_1_2[ind_1_2])
                ind_1_2 += 1
            else:
                try:
                    char_3 = int(s_with_dollar[sa_3[ind_3]])
                    char_1_2 = int(s_with_dollar[sa_1_2[ind_1_2]])
                except Exception as e:
                    char_3 = s_with_dollar[sa_3[ind_3]]
                    char_1_2 = s_with_dollar[sa_1_2[ind_1_2]]

                if char_3 < char_1_2:
                    merged.append(sa_3[ind_3])
                    ind_3 += 1
                elif char_3 > char_1_2:
                    merged.append(sa_1_2[ind_1_2])
                    ind_1_2 += 1
                elif sa_1_2[ind_1_2] % 3 == 1:
                    try:
                        rank_3 = sa_1_2_ranks[sa_3[ind_3] + 1]
                    except:
                        rank_3 = 0
                    # if ind_1_2 >= len(sa_1_2) -2:
                    #   rank_1_2 = 0

                    try:
                        rank_1_2 = sa_1_2_ranks[sa_1_2[ind_1_2] + 1]
                    except:
                        rank_1_2 = 0
                    if rank_3 < rank_1_2:
                        merged.append(sa_3[ind_3])
                        ind_3 += 1
                    else:
                        merged.append(sa_1_2[ind_1_2])
                        ind_1_2 += 1
                elif sa_1_2[ind_1_2] % 3 == 2:
                    try:
                        second_char_3 = int(s_with_dollar[sa_3[ind_3] + 1])
                        second_char_1_2 = int(s_with_dollar[sa_1_2[ind_1_2] + 1])
                    except:
                        second_char_3 = s_with_dollar[sa_3[ind_3] + 1]
                        second_char_1_2 = s_with_dollar[sa_1_2[ind_1_2] + 1]

                    if second_char_3 < second_char_1_2:
                        merged.append(sa_3[ind_3])
                        ind_3 += 1
                    elif second_char_3 > second_char_1_2:
                        merged.append(sa_1_2[ind_1_2])
                        ind_1_2 += 1
                    else:
                        try:
                            rank_3 = sa_1_2_ranks[sa_3[ind_3] + 2]
                        except:
                            rank_3 = 0
                        try:
                            rank_1_2 = sa_1_2_ranks[sa_1_2[ind_1_2] + 2]
                        except Exception as e:
                            # print(f'Error: {e}')
                            rank_1_2 = 0
                        if rank_3 < rank_1_2:
                            merged.append(sa_3[ind_3])
                            ind_3 += 1
                        else:
                            merged.append(sa_1_2[ind_1_2])
                            ind_1_2 += 1
                else:
                    raise Exception('Is this possible?')

            if ind_1_2 >= len(sa_1_2) and ind_3 >= len(sa_3):
                done = True

        return merged

    def do_it_naive(self, max_sort_len=None):
        suffixes = [i for i in range(len(self._s))]
        if max_sort_len is not None:
            suffixes.sort(key=lambda x: self._s[x:][:max_sort_len])  # + '$' + self._s[:x])
        else:
            suffixes.sort(key=lambda x: self._s[x:])  # + '$' + self._s[:x])
        return suffixes

    def do_it(self):
        if len(self._s) == 1:
            return [0]
        if len(self._s) == 0:
            return []

        triplet_1_2_names, triplet_1_2_inds, all_unique = self.construct_1_2()
        u, u_inds = self.make_u(triplet_1_2_names)
        if all_unique:
            sa_12 = self.reconstruct_from_inds(u, u_inds, triplet_1_2_inds)
            ranks_3 = self.get_ranks_3(self._s, sa_12)
            merged = self.merge(self._s + '$$$', ranks_3, sa_12)
            return merged
        else:
            if len(u) < 4:
                u_sorted = self.small_sort(u)
                sa_12 = self.reconstruct_from_sorted_u(u_sorted, u_inds, triplet_1_2_names, triplet_1_2_inds)
                ranks_3 = self.get_ranks_3(self._s, sa_12)
                merged = self.merge(self._s + '$$$', ranks_3, sa_12)
                return merged
            else:
                sa_u = SuffixArray(u)
                sa_u_merged = sa_u.do_it()
                sa_12 = []
                for i, el in enumerate(sa_u_merged):
                    if el >= len(u_inds):
                        pass  # $$$?
                    elif u_inds[el] == -1:  # hash sign from u
                        pass
                    else:
                        sa_12.append(u_inds[el])
                ranks_3 = self.get_ranks_3(self._s, sa_12)
                merged = self.merge(self._s + '$$$', ranks_3, sa_12)
                return merged

    def __init__(self, s):
        self._s = s
        # self.mod_1, self.mod_2, self.mod_3 = self.construct()


class SuffixArrayRankList(SuffixArray):
    """ Suffix Array using a list of integer ranks instead of a string

    eg: [97, 100, 97] instead of 'ada'
    Note: Can still be initialised with string (will automatically be converted to rank list).
    Used to construct suffix array in linear time with recursive dc3 algoirthm

    """

    def __init__(self, s):
        super().__init__(s)
        if type(s) == str:
            self._sr = [ord(ch) for ch in s]
        else:
            self._sr = s

    def get_item_at_ind(self, ind):
        return self._sr[ind]

    def get_suffix_at_ind(self, ind):
        return self._sr[ind:]

    # noinspection PyMethodMayBeStatic
    def get_triplet_at_ind(self, s_with_sentinel, ind):
        return tuple(s_with_sentinel[ind:ind+3])

    def assign_lex_names(self, triplets):

        curr_name = 0
        prev_val = '***'
        triplet_names = {}
        triplet_inds = {}
        all_unique = True

        for i, tup in enumerate(triplets):
            if tup[0] == prev_val:
                all_unique = False
            else:
                curr_name += 1
            triplet_names[tup] = curr_name
            triplet_inds[i] = tup
            prev_val = tup[0]

        return triplet_names, triplet_inds, all_unique

    def get_ranks_3(self, u, sa_12):
        ranks_1_2 = {el: i + 1 for i, el in enumerate(sa_12)}
        ranks_3 = [((u[i], 0), i) if i >= len(u) - 1 else ((u[i], ranks_1_2[i + 1]), i) for i in
                   range(0, len(u), 3)]
        ranks_3 = RadixSortSpecialList(ranks_3, all_same_len=True).sort()
        ranks_3_sa = [x[1] for x in ranks_3]
        return ranks_3_sa

    def make_u(self, triplets_with_names):

        u = []
        mod_1_inds = [i for i in range(1, len(self._sr), 3)]
        mod_2_inds = [i for i in range(2, len(self._sr), 3)]
        s_with_sentinel = self._sr + [0, 0]
        for ind in mod_1_inds:
            u.append(triplets_with_names[(self.get_triplet_at_ind(s_with_sentinel, ind), ind)])
        u.append(-1)
        for ind in mod_2_inds:
            u.append((triplets_with_names[(self.get_triplet_at_ind(s_with_sentinel, ind), ind)]))

        return u, mod_1_inds + [-1] + mod_2_inds

    def construct_1_2(self):

        s_with_sentinel = (self._sr + [0]) if len(self._sr) % 3 == 1 else (self._sr + [0, 0])
        # s_with_sentinel = (self._s + '$') if len(self._s) % 3 == 1 else (self._s + '$$')  # + '$$$'
        len_s = len(self._sr)
        triplets_1_2 = []
        for i in range(len_s):  # (len_s - 2):
            if i % 3 == 0:
                pass
            else:
                triplet = self.get_triplet_at_ind(s_with_sentinel, i)  # '*'.join(s_with_sentinel[i:i + 3])
                triplets_1_2.append((triplet, i))
        triplets_1_2 = RadixSortSpecialList(triplets_1_2, all_same_len=True).sort()

        triplet_1_2_names, triplet_1_2_inds, all_unique = self.assign_lex_names(triplets_1_2)

        return triplet_1_2_names, triplet_1_2_inds, all_unique

    def do_it(self):
        if len(self._sr) == 1:
            return [0]
        if len(self._sr) == 0:
            return []

        triplet_1_2_names, triplet_1_2_inds, all_unique = self.construct_1_2()
        u, u_inds = self.make_u(triplet_1_2_names)
        if all_unique:
            sa_12 = self.reconstruct_from_inds(u, u_inds, triplet_1_2_inds)
            ranks_3 = self.get_ranks_3(self._sr, sa_12)
            merged = self.merge(self._sr + [0, 0, 0], ranks_3, sa_12)
            return merged
        else:
            if len(u) < 4:
                u_sorted = self.small_sort(u)
                sa_12 = self.reconstruct_from_sorted_u(u_sorted, u_inds, triplet_1_2_names, triplet_1_2_inds)
                ranks_3 = self.get_ranks_3(self._sr, sa_12)
                merged = self.merge(self._sr + [0, 0, 0], ranks_3, sa_12)
                return merged
            else:
                sa_u = SuffixArrayRankList(u)
                sa_u_merged = sa_u.do_it()
                sa_12 = []
                for i, el in enumerate(sa_u_merged):
                    if el >= len(u_inds):
                        pass  # $$$?
                    elif u_inds[el] == -1:  # hash sign from u
                        pass
                    else:
                        sa_12.append(u_inds[el])
                ranks_3 = self.get_ranks_3(self._sr, sa_12)
                merged = self.merge(self._sr + [0, 0, 0], ranks_3, sa_12)
                return merged


class CountingSort:
    def __init__(self, items):
        self._items = items
        self.init_content_type()
        self.init_dict()

    def init_content_type(self):
        if type(self._items) == str:
            self.content_type = 's'
        else:
            self.content_type = 'l'

    def init_dict(self):
        self._counts = defaultdict(int)

    def sort(self):

        for item in self._items:
            self._counts[item] += 1

        sorted_keys = sorted(self._counts.keys())

        if self.content_type == 's':
            out_items = ''
            for key in sorted_keys:
                for i in range(self._counts[key]):
                    out_items += key
        else:
            out_items = []
            for key in sorted_keys:
                for i in range(self._counts[key]):
                    out_items.append(key)

        return out_items


class CountingSortExp(CountingSort):

    def __init__(self, items, exp=0, also_return_min_max=False, add_orig_ind=False):
        super().__init__(items)
        if add_orig_ind:
            self._items = [(item, i) for i, item in enumerate(self._items)]
            self.item_type = 't'
        self.exp = exp
        self.also_return_min_max = also_return_min_max

    def init_content_type(self):
        self.item_type = None
        if type(self._items) == list:
            self.content_type = 'l'
            if type(self._items[0]) == tuple:
                self.item_type = 't'  # tuple, also contains orig indexes
            else:
                self.item_type = 'n'  # normal
        else:
            raise Exception(f'Error - Must be list. Type provided is: {type(self._items)}')

    def init_dict(self):
        self._counts = defaultdict(list)  # Use list in order to preserve orig posns

    def _sort_string(self):

        min_num = math.inf
        max_num = math.inf * -1

        for item in self._items:
            if self.item_type == 't':
                # dig = item[0][len(item[0]) - 1 - self.exp]
                if len(item[0]) < self.exp + 1:
                    dig = ' '
                else:
                    dig = item[0][self.exp]
            else:
                # dig = item[len(item) - 1 - self.exp]
                if len(item) < self.exp + 1:
                    dig = ' '
                else:
                    dig = item[self.exp]
            self._counts[dig].append(item)
        return self._counts, min_num, max_num

    def _sort_numeric(self):

        min_num = math.inf
        max_num = math.inf * -1

        for item in self._items:
            if self.item_type == 't':
                dig = abs(item[0]) // pow(10, self.exp)
                ind = dig % 10 if item[0] >= 0 else (-1 * (dig % 10))
#                if item[0] < 0:
#                    ind = (dig * -1) % 10
#                    ind *= -1
#                else:
#                    ind = dig % 10

            else:
                dig = abs(item) // pow(10, self.exp)
                ind = dig % 10 if item >= 0 else (-1 * (dig % 10))
                # if item < 0:
                #     ind = (dig * -1) % 10
                #     ind *= -1
                # else:
                #     ind = dig % 10

            if self.also_return_min_max:
                if self.item_type == 't':
                    if item[0] < min_num:
                        min_num = item[0]
                    elif item[0] > max_num:
                        max_num = item[0]
                else:
                    if item < min_num:
                        min_num = item
                    elif item > max_num:
                        max_num = item
            self._counts[ind].append(item)
        return self._counts, min_num, max_num

    def sort(self):

        if self.item_type == 't':
            if type(self._items[0][0]) == str:
                item_content = 's'  # string
            else:
                item_content = 'n'  # number
        else:
            if type(self._items[0]) == str:
                item_content = 's'  # string
            else:
                item_content = 'n'  # number

        if item_content == 's':
            self._counts, min_num, max_num = self._sort_string()
        else:
            self._counts, min_num, max_num = self._sort_numeric()

        sorted_keys = sorted(self._counts.keys())

        out_items = []
        for key in sorted_keys:
            for item in self._counts[key]:
                out_items.append(item)

        if self.also_return_min_max:
            return out_items, min_num, max_num
        else:
            return out_items


class CountingSortExpSpecialList(CountingSortExp):

    def __init__(self, items, exp=0, also_return_min_max=False, add_orig_ind=False, assume_fixed_len=False):
        super().__init__(items, exp=exp, also_return_min_max=also_return_min_max, add_orig_ind=add_orig_ind)
        self.assume_fixed_len = assume_fixed_len
        self.fixed_len = None

        if len(items) == 0:
            pass
        else:
            if type(items) == list:
                pass
            else:
                raise Exception('Error: Must be list')
            if type(items[0]) == tuple:
                if type(items[0][0]) == tuple:
                    if self.assume_fixed_len:
                        self.fixed_len = len(items[0][0])
                elif type(items[[0][0] == int]):
                    if self.assume_fixed_len:
                        self.fixed_len = len(items[0])
                else:
                    raise Exception('Error: Must contain tuples of integer ranks')
            else:
                raise Exception('Error: Must contain tuples of integer ranks')

    def get_len(self, item):
        if self.fixed_len is not None:
            return self.fixed_len

        return len(item.split('*'))

    # noinspection PyMethodMayBeStatic
    def get_item_ind(self, item, ind):
        return int(item.split('*')[ind])

    def sort_rank_strings(self):
        min_num = math.inf
        max_num = math.inf * -1

        for item in self._items:
            if self.item_type == 't':
                # dig = item[0][len(item[0]) - 1 - self.exp]
                if self.get_len(item[0]) < self.exp + 1:
                    dig = 0
                else:
                    dig = self.get_item_ind(item[0], self.exp)
            else:
                # dig = item[len(item) - 1 - self.exp]
                if len(item) < self.exp + 1:
                    dig = ' '
                else:
                    dig = item[self.exp]
            self._counts[dig].append(item)
        return self._counts, min_num, max_num

    def sort_rank_lists(self):
        min_num = math.inf
        max_num = math.inf * -1

        for item in self._items:
            if self.item_type == 't':
                # dig = item[0][len(item[0]) - 1 - self.exp]
                if self.fixed_len is not None:
                    item_len = self.fixed_len
                else:
                    if type(item[0]) == tuple:
                        item_len = len(item[0])
                    else:
                        item_len = len(item)
                if item_len < self.exp + 1:
                    dig = 0
                else:
                    if type(item[0]) == tuple:
                        dig = item[0][self.exp]
                    else:
                        dig = item[self.exp]
            else:
                # dig = item[len(item) - 1 - self.exp]
                item_len = self.fixed_len if self.fixed_len is not None else len(item)
                if item_len < self.exp + 1:
                    dig = 0
                else:
                    dig = item[self.exp]
            self._counts[dig].append(item)
        return self._counts, min_num, max_num

    def sort(self):

        item_content = 'r'  # rank strings

        self._counts, min_num, max_num = self.sort_rank_lists()

        sorted_keys = sorted(self._counts.keys())

        out_items = []
        for key in sorted_keys:
            for item in self._counts[key]:
                out_items.append(item)

        if self.also_return_min_max:
            return out_items, min_num, max_num
        else:
            return out_items


class RadixSort:
    def __init__(self, items, add_orig_ind=False, all_same_len=False):
        self._items = items
        self.add_orig_ind = add_orig_ind
        self.all_same_len = all_same_len  # Assume all items same len. If not, need to pad

        if type(self._items[0]) == tuple:
            if type(self._items[0][0] == str):
                self.item_content = 's'
            else:
                self.item_content = 'n'
        else:
            if type(self._items[0]) == str:
                self.item_content = 's'
            else:
                self.item_content = 'n'

    def sort(self):
        if self.item_content == 's':
            if self.all_same_len:
                if type(self._items[0]) == tuple:
                    max_len = len(self._items[0][0])
                else:
                    max_len = len(self._items[0])
            else:
                if type(self._items[0]) == tuple:
                    lens = [len(s[0]) for s in self._items]
                else:
                    lens = [len(s) for s in self._items]
                max_len = max(lens)

            # items = CountingSortExp(self._items, exp=0, add_orig_ind=self.add_orig_ind).sort()
            items = CountingSortExp(self._items, exp=max_len-1, add_orig_ind=self.add_orig_ind).sort()

            # for exp in range(1,  max_len): # len(items[0][0])):
            #   items = CountingSortExp(items, exp=exp).sort()

            for exp in range(max_len - 2, -1, -1):  # len(items[0][0])):
                items = CountingSortExp(items, exp=exp).sort()

        else:

            items, min_num, max_num = CountingSortExp(self._items, exp=0, also_return_min_max=True,
                                                      add_orig_ind=self.add_orig_ind).sort()
            if max_num >= 0 and min_num >= 0:
                log = int(math.log(max_num, 10))
            elif max_num < 0 and min_num < 0:
                log = int(math.log(min_num*-1, 10))
            else:
                max_log = int(math.log(max_num, 10))
                min_log = int(math.log(min_num*-1, 10))
                log = max(min_log, max_log)

            for exp in range(1, log+1):
                items = CountingSortExp(items, exp=exp).sort()

        return items


class RadixSortSpecialList(RadixSort):

    def __init__(self, items, add_orig_ind=False, all_same_len=False):
        super().__init__(items, add_orig_ind=add_orig_ind, all_same_len=all_same_len)

        self.item_content = 'r'  # ranks

        if len(items) == 0:
            pass
        else:
            if type(items) == list:
                pass
            else:
                raise Exception('Error: Must be list')
            if type(items[0]) == tuple:
                if type(items[0][0]) == tuple:
                    pass  # Each item contains tuple of (ranks tuple, ind)
                elif type(items[0][0]) == int:
                    pass  # Each item contains (ranks tuple)
                else:
                    raise Exception('Error: Tuples must contain tuple of ranks as 1st element')
            else:
                raise Exception('Error: Must be list of rank tuples')

    # noinspection PyMethodMayBeStatic
    def get_item_len(self, item):
        return len(item.split('*'))

    def sort(self):
        if self.item_content != 'r':
            return super.sort()  # Just ordinary radix for list of strings or list of numbers

        if self.all_same_len:
            if type(self._items[0][0]) == tuple:
                max_len = len(self._items[0][0])
            else:
                max_len = len(self._items[0])
        else:
            if type(self._items[0][0]) == tuple:
                lens = [len(s[0]) for s in self._items]
            else:
                lens = [len(s) for s in self._items]
            max_len = max(lens)

        # items = CountingSortExp(self._items, exp=0, add_orig_ind=self.add_orig_ind).sort()
        items = CountingSortExpSpecialList(self._items, exp=max_len-1, add_orig_ind=self.add_orig_ind,
                                           assume_fixed_len=self.all_same_len).sort()

        # for exp in range(1,  max_len): # len(items[0][0])):
        #   items = CountingSortExp(items, exp=exp).sort()

        for exp in range(max_len - 2, -1, -1):  # len(items[0][0])):
            items = CountingSortExpSpecialList(items, exp=exp, assume_fixed_len=self.all_same_len).sort()

        return items


def gen_rand_string(length, alphabet='ACGT'):
    str_list = []
    for i in range(length):
        r = random.randint(0, len(alphabet)-1)
        str_list.append(alphabet[random.randint(0, len(alphabet)-1)])
    return ''.join(str_list)


if __name__ == '__main__':

    #    s_in = 'panamabananas'

    str_len = None  # 100000 #None #500000
    sequences = read_fasta('../test_data/large/GCF_000005845.2_ASM584v2_genomic.fna')
    if str_len is None:
        str_len = len(sequences[0][1])
        print(f'Using full string: {str_len}')
    print(f'ecoli suffix array: {str_len} ')
    s_in = sequences[0][1][:str_len]

    # s_in = gen_rand_string(100000)
    # print(s_in)
    start = time.time()
    sa = SuffixArrayRankList(s_in)
    res_l = sa.do_it()
    end = time.time()
    print(f'RankList: {end - start}')
    exit()

    start = time.time()
    res_naive = sa.do_it_naive()
    for ii, el in enumerate(res_l):
        if el != res_naive[ii]:
            print(f'error at pos {ii}')
            break
    end = time.time()
    print(f'Naive: {end - start}')

    if False:

        str_len = None  # 100000 #None #500000
        sequences = read_fasta('../test_data/large/GCF_000005845.2_ASM584v2_genomic.fna')
        if str_len is None:
            str_len = len(sequences[0][1])
            print(f'Using full string: {str_len}')
        print(f'ecoli suffix array: {str_len} ')
        s_ecoli = sequences[0][1][:str_len]
        start = time.time()
        sa = SuffixArrayRankList(s_ecoli)
        res_l = sa.do_it()
        end = time.time()
        print(f'  RankList: {end - start}')

        print(res_l[:10], s_ecoli[res_l[0]:res_l[0]+20], s_ecoli[res_l[1]:res_l[1] + 20],
              s_ecoli[res_l[-1]:res_l[-1] + 20])
        # start = time.time()
        # res2 = sa.do_it_naive(max_sort_len=None)
        # end = time.time()
        # print(end - start)

        start = time.time()
        res_naive = sa.do_it_naive(max_sort_len=None if str_len < 200 else 200)
        for ii, el in enumerate(res_l):
            if el != res_naive[ii]:
                print(f'error at pos {ii}')
                break
        end = time.time()
        print(f'  Naive: {end - start}')

        exit()
    s_in = 'panamabananas'
    bw = BWT(s_in, verbose=True, use_linear_sa_construction=True)
    t_in = 'ana'
    matches = bw.pattern_match(t_in, force_reconstruct=True)
    print(matches)

    s_in = gen_rand_string(5000)
    bw = BWT(s_in, verbose=1, use_linear_sa_construction=False, max_sort_len_for_naive=10000)
    t_in = 'ACGT'
    matches = bw.pattern_match(t_in, force_reconstruct=True)
    print(len(matches))
    bw = BWT(s_in, verbose=1, use_linear_sa_construction=True)
    t_in = 'ACGT'
    matches = bw.pattern_match(t_in, force_reconstruct=True)
    print(len(matches))

    xx = 1

    if False:
        # Construct bwt for full ecoli
        str_len = None
        sequences = read_fasta('../test_data/GCF_000005845.2_ASM584v2_genomic.fna')
        if str_len is None:
            str_len = len(sequences[0][1])
        s_in = sequences[0][1][:str_len]
        t_in = 'GATTACA'
        bw = BWT(s_in, verbose=1)
        res = bw.pattern_match(t_in)
        print(len(res))
        print(res[:5])
        list_to_file(f'../test_data/ecoli_{str_len}.sa', bw.sa)
        string_to_file(f'../test_data/ecoli_{str_len}.bw', bw.bw)
        sa_recovered = file_to_list(f'../test_data/ecoli_{str_len}.sa')
        bw_recovered = read_txt_file(f'../test_data/ecoli_{str_len}.bw')
        bw = BWT.from_bw(bw_recovered, sa_recovered)
        res = bw.pattern_match(t_in)
        print(len(res))
        print(res[:5])

        exit()

    str_len = 4641652
    sa_recovered = file_to_list(f'../test_data/ecoli_{str_len}.sa')
    bw_recovered = read_txt_file(f'../test_data/ecoli_{str_len}.bw')
    bw = BWT.from_bw(bw_recovered, sa_recovered, verbose=1)
    t_in = 'GATTACA'
    res = bw.pattern_match(t_in)
    print(len(res))
    print(res[:5])
    sequences = read_fasta('../test_data/GCF_000005845.2_ASM584v2_genomic.fna')
    s_orig = sequences[0][1][:str_len]

    if False:
        str_len = None  # 100000
        sequences = read_fasta('../test_data/GCF_000005845.2_ASM584v2_genomic.fna')
        if str_len is None:
            str_len = len(sequences[0][1])
            print(f'Using full string: {str_len}')
        s_ecoli = sequences[0][1][:str_len]
        sa = SuffixArrayRankString(s_ecoli)
        start = time.time()
        res = sa.do_it()
        end = time.time()
        print(end - start)
        print(res[:10], s_ecoli[res[0]:res[0]+20], s_ecoli[res[1]:res[1] + 20], s_ecoli[res[-1]:res[-1] + 20])
        # start = time.time()
        # res2 = sa.do_it_naive(max_sort_len=None)
        # end = time.time()
        # print(end - start)

        list_to_file(f'../test_data/ecoli_{str_len}.sa', res)

        exit()

    rs = RadixSortSpecial(['11*12*12*11*12', '14*7', '12*0', '12*9*3', '12', '11*12', '12*9', '12*-8', '12*9', '12*0',
                           '11'], add_orig_ind=True)
    res2 = rs.sort()

    small_sample_s = 'AAAAGTCGCTGCAGCGT'
    sample_s = 'AAAAGTCGCTGCAGCGTCGAGAGAGCAAAAAAATTAGGCGATGCGGAGC'

    rs = RadixSortSpecial(['12', '11*12', '12*11*12', '12*12*11*12', '11*12*12*11*12'], add_orig_ind=True)
    res3 = rs.sort()

    sa = SuffixArray('banana')
    res4 = sa.do_it()
    sa = SuffixArrayRankString('12*11*13*11*13*11', rank_string_input=True)
    res5 = sa.do_it()

    sa = SuffixArray('panamabananas')
    res6 = sa.do_it()
    sa = SuffixArrayRankString('15*11*14*11*13*11*12*11*14*11*14*11*16', rank_string_input=True)
    res7 = sa.do_it()
    sa = SuffixArrayRankString('112*97*110*97*109*97*98*97*110*97*110*97*115', rank_string_input=True)
    res8 = sa.do_it()
    sa = SuffixArrayRankString('panamabananas')
    res9 = sa.do_it()

    s_in = 'GCAGTGA'
    sa = SuffixArrayRankString(s_in)
    res10 = sa.do_it()
    res11 = sa.do_it_naive()
    for ii, ch_ch in enumerate(res):
        if ch_ch != res2[ii]:
            print(f'aaargh: {ii} {s_in}')
            break

    for ii in range(10):
        s_len = random.randint(1, 100)
        s_in = gen_rand_string(s_len)
        print(s_in)
        sa = SuffixArrayRankString(s_in)
        res12 = sa.do_it()
        res13 = sa.do_it_naive()
        for iii, ch_ch in enumerate(res):
            if ch_ch != res13[iii]:
                print(f'aaargh: {iii} {s_in}')
                break

    big_sample_s = 'AATACCGCTGAGGTAGCGACTTAGACATAAAAGAGGGGCAATGCCAACCGGTATCTTTATCCCCGACACGCTCGGTCCGATAGGGCGCTACTCAATCAAGGCTGGCAGCCAGTGCAAGTGATCTTTCAGTGCTAGCGTAATGGTTTGCGCCCGAGTAAGTCTGAAGCCCGCTCCGCTATCGATATACCGCGTCGATGGGGTTACTTAAAACTGTAGCTAGGATATTCCGTAGGGGCCCGGCTAGGGCGCATGAATCAATCTCACCGATCCAGGAGAAGAAGACTGGTCCCCAGCCCGATC'

    very_big_sample_s = gen_rand_string(500000)
    sa = SuffixArrayRankString(very_big_sample_s)
    start = time.time()
    res14 = sa.do_it()
    end = time.time()
    print(end - start)
    # start = time.time()
    # res2 = sa.do_it_naive()
    # end = time.time()
    # print(end - start)

    Z = pattern_match_zbox("aaaaaa", t_in=None, also_return_zbox=True)  # X 5 4 3 2 1
    print(Z)
    Z = pattern_match_zbox("aabaacd", t_in=None, also_return_zbox=True)   # X 1 0 2 1 0 0
    print(Z)
    Z = pattern_match_zbox("abababab", t_in=None, also_return_zbox=True)  # X 0 6 0 4 0 2 0
    print(Z)
    Z = pattern_match_zbox("ACACACGTACACG", t_in=None, also_return_zbox=True)  # X 0 4 0 2 0 0 0 4 0 2 0 0
    print(Z)
    Z = pattern_match_zbox(sample_s, 'AAG', also_return_zboxTrue)
    print(Z)
    repeat_s = 'AAAGTAAAGGCAAAGAAAGCCAAAAG' * 1000
    Z = pattern_match_naive(repeat_s, 'AAAG')
    print(Z)
    Z = pattern_match_zbox(repeat_s, 'AAAG', also_return_zbox=True)
    print(Z)
    # worst case for naive:
    s_in = 'a' * 100000
    t_in = 'a' * 100
    Z = pattern_match_zbox(s_in, t_in, also_return_zbox=True)
    print(Z)
    Z = pattern_match_naive(s_in, t_in)
    print(Z)

    s_in = 'panamabananas'
    bw = BWT(s_in)
    orig = bw.reconstruct_naive(also_return_suffix_array=True)
    print(orig)
    res15 = bw.pattern_match('ana')
    bw.save('../test_data/panamabananas.bw', also_save_suffix_array=True)

    bw = BWT.from_bw('enwvpeoseu$llt')
    orig = bw.reconstruct_naive(also_return_suffix_array=True)
    res16 = bw.pattern_match('one')
    res17 = bw.pattern_match('plus', force_reconstruct=True)
    res18 = bw.pattern_match('one')

    ecoli = True
    if ecoli:
        sequences = read_fasta('../test_data/GCF_000005845.2_ASM584v2_genomic.fna')
        s_in = sequences[0][1][:500000]
        t_in = s_in[100000:101000]  # 'GATGATGAT'
        # t_in = 'GATTACA'
        print(len(s_in), len(t_in), len(s_in) * len(t_in), len(s_in) + len(t_in))
        matches = pattern_match_naive(s_in, t_in)
        print(len(matches))
        matches = pattern_match_zbox(s_in, t_in)
        print(len(matches))
        matches = pattern_match_zbox(s_in, t_in)
        print(len(matches))
        bw_s = read_txt_file('../test_data/ecoli_500000.bw')
        bw_sa = file_to_list('../test_data/ecoli_500000.sa')
        bw = BWT.from_bw(bw_s, sa=bw_sa)
        # bw = BWT(s_in)
        matches = bw.pattern_match(t_in, force_reconstruct=True)
        print(len(matches))
        matches = bw.pattern_match(t_in)
        print(len(matches))

    # string_to_file('../test_data/ecoli_500000.bw', bw.bw)
    # xx = 1/0

    # sequences = read_fasta('../test_data/pseudo_dna.fa')
    # s_in = sequences[0][1]
    s_in = read_txt_file('../test_data/pseudo_dna_seq.bw')
    t_in = s_in[100:107]
    bw = BWT(s_in)
    # bw = BWT.from_bw(s_in)
    # s_prime, s_prime_inds = bw.reconstruct_naive()
    print(s_in[:20], bw.bw[:20])
    res19 = bw.pattern_match(t_in)
    res20 = bw.pattern_match(t_in, force_reconstruct=True)

    comp = bw.compress_bw()

    # string_to_file('../test_data/pseudo_dna_seq.bw', bw.bw)

    s_in = sample_s  # 'panamabananas'
    bw = BWT(s_in)
    s_prime = bw.reconstruct_naive()
    print(s_in, bw.bw, s_prime)

    res21 = bw.pattern_match('GAG')  # ('ana')
    # matches = pattern_match_naive(sample_s, 'AAA')
    # print(matches)

    haystack = gen_rand_string(10000, alphabet='deln')
    haystack = haystack[:100] + 'needle' + haystack[100:]
    bw = BWT(haystack)
    res22 = bw.pattern_match('needle')

    # s_in = gen_rand_string(1000) #1000000
    # s_in = s_in * 10 # 1000
    sequences = read_fasta('../test_data/pseudo_dna.fa')
    s_in = sequences[0][1]
    bw_s = read_txt_file('../test_data/pseudo_dna_seq.bw')
    bw = BWT.from_bw(bw_s)
    t_in = 'GATTACA'
    # bw = BWT(s_in, also_return_match_inds=False)

    # t_in = s_in[30:100] #s_in[1000000:2000000]
    print(len(s_in), len(t_in), len(s_in) * len(t_in), len(s_in) + len(t_in))
    matches = pattern_match_naive(s_in, t_in)
    print(len(matches))
    matches = pattern_match_zbox(s_in, t_in)
    print(len(matches))
    bw2 = BWT(s_in)
    matches = bw2.pattern_match(t_in)
    print(len(matches))

    Z = pattern_match_zbox(sample_s, 'AAG', print_z=True)
    print(Z)
    Z = pattern_match_zbox("aaaaaa", None, print_z=True)  # X 5 4 3 2 1
    print(Z)
    Z = pattern_match_zbox("aabaacd", None, print_z=True)   # X 1 0 2 1 0 0
    print(Z)
    Z = pattern_match_zbox("abababab", None, print_z=True)  # X 0 6 0 4 0 2 0
    print(Z)
    Z = pattern_match_zbox("ACACACGTACACG", None, print_zTrue)  # X 0 4 0 2 0 0 0 4 0 2 0 0
    print(Z)
