import numpy as np
import json


class PopDistGen:
    """Class representing Population Distribution of a single generation"""

    def __init__(self, fa, pop, gen_num, genotype_fitnesses=[1.0, 1.0, 1.0], F=0, verbose=0, prev_gen_num=None):
        """

        Args:
            fa:
            pop:
            gen_num:
            genotype_fitnesses:
            F:
            verbose:
            prev_gen_num:
        """
        self.in_fa = fa
        self.pop = pop
        self.genotype_fitnesses = genotype_fitnesses
        self.F = F
        self.gen_num = gen_num
        self.prev_gen_num = prev_gen_num
        self.verbose = verbose

        self.random_mating_genotypes = None
        # self.survived_genotypes = None

        if prev_gen_num is None:
            self.out_fa = fa
        else:
            self.perform_matings()

    @property
    def in_fA(self):
        return 1 - self.in_fa

    @property
    def out_fA(self):
        return 1 - self.out_fa

    # noinspection PyMethodMayBeStatic
    def calc_adjust_factor(self, raw_freqs):
        """ calculate adjustment to maintain population size"""
        return 1.0 / sum(raw_freqs)

    @property
    def survived_genotypes(self):

        survived_AA = round(self.survived_genotype_freqs[0] * self.pop)
        survived_Aa = round(self.survived_genotype_freqs[1] * self.pop)
        survived_aa = self.pop - survived_AA - survived_Aa
        return [survived_AA, survived_Aa, survived_aa]

    def survival(self):
        raw_survived_genotype_freqs = [self.genotype_fitnesses[i] * genotype_freq
                                       for i, genotype_freq in enumerate(self.random_mating_genotype_freqs)]
        adjust_factor = self.calc_adjust_factor(raw_survived_genotype_freqs)

        adjusted_genotype_freqs = [genotype_freq * adjust_factor for genotype_freq in raw_survived_genotype_freqs]

        if self.verbose > 0:
            print('mating ', self.gen_num, ' rand: ', self.random_mating_genotype_freqs, 'raw survived: ',
                  raw_survived_genotype_freqs, ' survived: ', adjusted_genotype_freqs)
        return adjusted_genotype_freqs

    def perform_matings(self):

        if self.pop is None:
            self.random_mating_genotype_freqs = PopDistGen.hw_genotype_freqs(self.in_fa, F=self.F)
        else:
            self.random_mating_genotype_counts = list(
                np.random.multinomial(self.pop, PopDistGen.hw_genotype_freqs(self.in_fa, F=self.F)))
            self.random_mating_genotype_counts = [int(rmgc) for rmgc in self.random_mating_genotype_counts]
            self.random_mating_genotype_freqs = [count / self.pop for count in self.random_mating_genotype_counts]

        self.survived_genotype_freqs = self.survival()
        self.out_fa = PopDistGen.allele_freqs_from_genotype_freqs(self.survived_genotype_freqs)[1]
        if self.verbose > 0:
            print(self)

    @staticmethod
    def hw_genotype_freqs(fa, F=0):

        fpq = F * fa * (1-fa)
        return (1-fa)**2 + fpq, 2*(fa*(1-fa) - fpq), fa**2 + fpq

    @staticmethod
    def allele_freqs_from_genotype_freqs(genotype_freqs):
        return genotype_freqs[0] + 0.5*genotype_freqs[1], genotype_freqs[2] + 0.5*genotype_freqs[1]

    def __str__(self):
        return 'PopDistGen: In fA/fa: %f %f Pop: %s, Gen: %i, Out fA/fa: %f %f'\
               % (self.in_fA, self.in_fa, self.pop, self.gen_num, self.out_fA, self.out_fa)

    @staticmethod
    def pop_dist_gen_from_prev_gen(prev_gen):
        return PopDistGen(prev_gen.out_fa, prev_gen.pop, prev_gen.gen_num + 1,
                          genotype_fitnesses=prev_gen.genotype_fitnesses, F=prev_gen.F,
                          verbose=prev_gen.verbose, prev_gen_num=prev_gen.gen_num)

    def to_json(self):
        # return '{ "init_fa": ' + self.init_fa + ', "pop": ' self.pop + '}'
        return vars(self)


class PopDist:
    """Class representing Population Distribution"""
    def __init__(self, fa, genotype_fitnesses=[1.0, 1.0, 1.0], pop=None, F=0, verbose=0):
        """ Constructor
         **Bold** statement of an *italic*

         List:
         - One
         - Two
         - Three

         > Note

         """
        self.init_fa = fa
        self.init_pop = pop
        self.init_genotype_fitnesses = genotype_fitnesses
        self.init_F = F
        self.verbose = verbose

        init_gen = PopDistGen(self.init_fa, self.init_pop, 0, genotype_fitnesses=self.init_genotype_fitnesses,
                              F=self.init_F, verbose=self.verbose)
        print(init_gen)

        self.gens = []
        self.gens.append(init_gen)

    @property
    def init_fA(self):
        return 1 - self.init_fa

    @property
    def num_gens(self):
        return len(self.gens) - 1

    @property
    def current_fA(self):
        return self.gens[self.num_gens].out_fA

    @property
    def current_fa(self):
        return self.gens[self.num_gens].out_fa

    def current_allele_freqs(self):
        return self.current_fA, self.current_fa

    def get_prev_gen(self, curr_gen):
        """
         Get previous generation

        Args:
           curr_gen: current generation

        Returns:
           previous generation

        """
        return self.gens[curr_gen-1]

    def sim_generations(self, num_generations):
        for i in range(1, num_generations+1):
            pdg = PopDistGen.pop_dist_gen_from_prev_gen(self.get_prev_gen(i))
            self.gens.append(pdg)
            # print('pdg: ',pdg)

    def __str__(self):
        return 'Init fA/fa: %f %f Curr fA/fa: %f %f' % (self.init_fA, self.init_fa, self.current_fA, self.current_fa)

    def get_plot_data(self, allele=1):
        # alleles: 0 - A, 1 - a, 2 - both
        # TODO  (both alleles not implemented yet)
        print('allele choice: ', allele)
        x_data = [i for i in range(len(self.gens))]
        y_data = [round(gen.out_fa if allele == 1 else gen.out_fA, 3) for gen in self.gens]
        return {'x_data': x_data, 'y_data': y_data}

#    @staticmethod
#    def Genotype_Freq_From_Allele_Freq(fA,fa):

    def _attr_dict(self):
        attr_dict = vars(self)
        attr_dict_out = attr_dict.copy()
        attr_dict_out['gens'] = [vars(gen) for gen in self.gens]
        return attr_dict_out

    def to_json(self):
        return json.dumps(self._attr_dict())

    @staticmethod
    def pop_dist_from_json(j):
        j_dict = json.loads(j)
        pop_dist = PopDist(j_dict['init_fa'], j_dict['init_genotype_fitnesses'], j_dict['init_pop'],
                           j_dict['init_F'], j_dict['verbose'])
        pop_dist.gens.pop()  # remove default init gen
        for j_gen in j_dict['gens']:
            pop_dist_gen = PopDistGen(j_gen['in_fa'], j_gen['pop'], j_gen['gen_num'],
                                      j_gen['genotype_fitnesses'], j_gen['F'], j_gen['verbose'], j_gen['prev_gen_num'])
            pop_dist_gen.out_fa = j_gen['out_fa']
            pop_dist.gens.append(pop_dist_gen)
        return pop_dist


if __name__ == '__main__':
    np.random.seed(42)
    n = 10
    p = 0.5
    nums = []
    for i in range(0, 100):
        nums.append(np.random.binomial(n, p))
        print(nums[i])

    pdhw = PopDist(0.2, pop=1000, F=0.4, verbose=1)
    print(pdhw)
    pdhw.sim_generations(1)
    print(PopDistGen.hw_genotype_freqs(0.2, F=0.4))

    pd = PopDist(0.5, genotype_fitnesses=[0.9, 1.0, 1.0], pop=100, F=0, verbose=1)
    print(pd)

    pd.sim_generations(3)

    j_out = pd.to_json()

    pd_inflated = PopDist.pop_dist_from_json(j_out)
    print('inf: ', pd_inflated)
    print(pd_inflated.gens[1])
