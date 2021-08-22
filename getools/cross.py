"""
Classes used to define organisms and perform organism crosses
"""

import random
import math
import itertools
import re

debug = 0
from getools.utils.gen_utils import SerialiserMixin
from scipy.stats import chisquare
import random

class Allele(SerialiserMixin):
    """
    Represents an individual allele (dominant or recessive - defaults to recessive)
    """

    Dominant = 'D'
    Recessive = 'R'

    def __init__(self,symbol,type='R'):
        """
        Constructor

        Args:
            symbol (str) : symbol used to represent allele (eg "A")
            type (str, optional) : Allele type (eg "R" - recessive or "D" - dominant). Default "R"
        """
        self.symbol = symbol
        self.type = type

    def __str__(self):
        return self.symbol

    def is_dom(self):
        return self.type == Allele.Dominant

    def is_rec(self):
        return self.type == Allele.Recessive

    @staticmethod
    def _from_attr_dict(attr_dict):
        obj = Allele(attr_dict['symbol'],attr_dict['type'])
        return obj

class AlleleSet(SerialiserMixin):
    """
    Set of alleles which are possible for a gene
    """
    def __init__(self,alleles):
        self.alleles= alleles

    @staticmethod
    def default_alleleset_from_symbol(symbol):
        alleles = []
        allele = Allele(symbol.upper(),type=Allele.Dominant)
        alleles.append(allele)
        allele = Allele(symbol.lower(),type=Allele.Recessive)
        alleles.append(allele)
        return AlleleSet(alleles)

    def dom_allele(self):
        for allele in self.alleles:
            if allele.is_dom():
                return allele
        return None

    def rec_allele(self):
        for allele in self.alleles:
            if allele.is_rec():
                return allele
        return None

    def __str__(self):
        return '-'.join([str(al) for al in self.alleles])



    @staticmethod
    def _from_attr_dict(attr_dict):
        alleles = [Allele._from_attr_dict(al) for al in attr_dict['alleles']]
        obj = AlleleSet(alleles)
        return obj

class Gene(SerialiserMixin):

    """
    Consists of an AlleleSet, a position, a name and an inheritance pattern
    """

    INH_PATTERN_DOMINANT= 'D'
    INH_PATTERN_RECESSIVE = 'R'

    def __init__(self,alleleset,position,name = '', inheritance_pattern=INH_PATTERN_RECESSIVE):
        self.alleleset = alleleset
        self.position = position
        self.name = name
        self.inheritance_pattern = inheritance_pattern
        #self.recessive = recessive

    def __str__(self):
        return 'Name: ' + self.name + ' Pos: ' +  str(self.position) + ' Alleles: ' + str(self.alleleset) + ' IP: ' + str(self.inheritance_pattern)

    def distance(self,other_gene):
        return other_gene.position - self.position

    @staticmethod
    def _from_attr_dict(attr_dict):
        alleleset = AlleleSet._from_attr_dict(attr_dict['alleleset'])
        obj = Gene(alleleset,attr_dict['position'],attr_dict['name'])
        return obj


class ChromosomeTemplate(SerialiserMixin):
    """ Template for a chromosome type

    Each actual chromosome belonging to an organism is linked to a specific ChromosomeTemplate to indicate its type
    (eg Automosomal / X / Y, size, genes associated with the chromosome etc)
    """

    #Class constants
    AUTOSOMAL = 'A'
    X = 'X'
    Y = 'Y'

    def __init__(self,name, size, genes_list=None, type=AUTOSOMAL):
        self.name = name
        self.size  = size
        self.type = type
        self.genes = self.add_genes(genes_list)

    @staticmethod
    def from_symbol_list(symbol_list, name, size=None, type=None):
        genes = []
        if type is None:
           type = ChromosomeTemplate.AUTOSOMAL
        default_position = ChromosomePair.UNLINKED_THRESHOLD
        if size is None:
           size = (len(symbol_list) + 1) * default_position

        for i, symbol_string in enumerate(symbol_list):
            split = symbol_string.split('-')
            symbol = split[0]
            inh_pattern = Gene.INH_PATTERN_RECESSIVE
            position = default_position * (i+1)
            if len(split) > 1:
               inh_pattern = split[1]
            if len(split) > 2:
               position = int(split[2])
            gene = Gene(AlleleSet.default_alleleset_from_symbol(symbol), position=position, inheritance_pattern=inh_pattern)
            genes.append(gene)
        return ChromosomeTemplate(name, size=size, genes_list=genes, type=type)



    def add_genes(self,genes_list):
         return sorted(genes_list,key=lambda x:x.position)

    def positions(self):
        return [str(gene.alleleset.alleles[1]) + '-' + str(gene.position) for gene in self.genes]

    def __str__(self):
        return 'ChromosomeTemplate ' + self.name +  '-Type: ' + self.type +  ' Size: ' + str(self.size) + ' Genes: ' + ','.join(['[Gene: ' + str(gene) + ']' for gene in self.genes])

    @staticmethod
    def _from_attr_dict(attr_dict):
        genes_list = [Gene._from_attr_dict(g) for g in attr_dict['genes']]
        obj = ChromosomeTemplate(attr_dict['name'],attr_dict['size'],genes_list)
        return obj

    def get_gene_from_symbol(self, symbol):
        for gene in self.genes:
            if symbol in str(gene.alleleset):
                return gene
        return None


    def generate_random_chromosome(self):
        actual_alleles = []
        for gene in self.genes:
            r = random.randint(0, len(gene.alleleset.alleles) - 1)
            actual_alleles.append(gene.alleleset.alleles[r])
        return Chromosome(self, actual_alleles)

    def generate_complement_chromosome(self, other_chromosome):
        actual_alleles = []
        for i,gene in enumerate(self.genes):
            other_allele = other_chromosome.alleles[i]
            other_allele_index = gene.alleleset.alleles.index(other_allele)
            allele_index = 0 if other_allele_index == 1 else 1
            actual_alleles.append(gene.alleleset.alleles[allele_index])
        return Chromosome(self, actual_alleles)


    def generate_hom_recessive_chromosome(self):
        actual_alleles = []
        for gene in self.genes:
            actual_alleles.append(gene.alleleset.rec_allele())  #(gene.alleles[-1])
        return Chromosome(self, actual_alleles)

    def generate_hom_dominant_chromosome(self):
        actual_alleles = []
        for gene in self.genes:
            actual_alleles.append(gene.alleleset.dom_allele())  #(gene.alleles[0])
        return Chromosome(self, actual_alleles)

    def generate_het_chromosome(self):
        actual_alleles = []
        for gene in self.genes:
            actual_alleles.append(gene.alleleset.alleles[-1])
        return Chromosome(self, actual_alleles)

    def generate_random_pair(self,ploidy):
        chrom_pair = []
        for i in range(ploidy):
            chromosome = self.generate_random_chromosome()
            chrom_pair.append(chromosome)
        return ChromosomePair(self,chrom_pair)


    def generate_hom_recessive_pair(self, ploidy):
        chrom_pair = []
        for i in range(ploidy):
            chromosome = self.generate_hom_recessive_chromosome()
            chrom_pair.append(chromosome)
        return ChromosomePair(self, chrom_pair)

    def generate_hom_dominant_pair(self, ploidy):
        chrom_pair = []
        for i in range(ploidy):
            chromosome = self.generate_hom_dominant_chromosome()
            chrom_pair.append(chromosome)
        return ChromosomePair(self, chrom_pair)

    def generate_het_pair(self, ploidy, rand_phase=False):
        chrom_pair = []
        if rand_phase:
            chromosome_1 = self.generate_random_chromosome()
            chromosome_2 = self.generate_complement_chromosome(chromosome_1)
        else:
            chromosome_1 = self.generate_hom_dominant_chromosome()
            chromosome_2 = self.generate_hom_recessive_chromosome()
        chrom_pair = [chromosome_1,chromosome_2]
        return ChromosomePair(self, chrom_pair)


class GenomeTemplate(SerialiserMixin):
    """ Template for a genome type

    Each actual organism's genome has a specific GenomeTemplate associated with it to indicate its type
    (eg ploidy, which chromosome templates it has, name)
    """
    def __init__(self,  ploidy=2, chromosome_templates=[],X_chromosome_template=None, Y_chromosome_template=None, name='Unnamed_Genome'):
        self.ploidy = ploidy
        self.chromosome_templates = chromosome_templates # autosomal chromosome_templates
        self.X_chromosome_template = X_chromosome_template
        self.Y_chromosome_template = Y_chromosome_template
        self.name = name

    def __str__(self):

        if self.X_chromosome_template is None:
            X_chrom = 'None'
        else:
            X_chrom = str(self.X_chromosome_template)


        if self.Y_chromosome_template is None:
            Y_chrom = 'None'
        else:
            Y_chrom = str(self.Y_chromosome_template)

        return 'GenomeTemplate: ' + self.name +  ':'.join(['[' + str(ct) + ']' for ct in self.chromosome_templates]) + ' X chrom: ' + X_chrom + ' Y chrom: ' + Y_chrom

    @staticmethod
    def _from_attr_dict(attr_dict):
        chromosome_templates  = [ChromosomeTemplate._from_attr_dict(ct) for ct in attr_dict['chromosome_templates']]
        obj = GenomeTemplate(ploidy=attr_dict['ploidy'],chromosome_templates=chromosome_templates, name=attr_dict['name'])
        return obj

    def positions(self):
        posns = []
        for chromosome_template in self.chromosome_templates:
            posns.extend(chromosome_template.positions())

        return posns

    def chromosome_names(self):
        return [chrom.name for chrom in self.chromosome_templates]

    def generate_hom_dominant_sex_pair(self, sex=None):
        if self.X_chromosome_template is None:
            return None

        chrom_1 = self.X_chromosome_template.generate_hom_dominant_chromosome() # at least 1 X
        if (self.Y_chromosome_template is None) or (sex == Genome.FEMALE):
            r = 0
        elif sex == Genome.MALE:
            r = 1
        else:
            r = random.randint(0,1)
        chrom_2_template = self.X_chromosome_template if r == 0 else self.Y_chromosome_template
        chrom_2 = chrom_2_template.generate_hom_dominant_chromosome()
        return ChromosomeSexPair([chrom_1,chrom_2], self.X_chromosome_template, chrom_2_template)

    def generate_hom_recessive_sex_pair(self, sex=None):
        if self.X_chromosome_template is None:
            return None

        chrom_1 = self.X_chromosome_template.generate_hom_recessive_chromosome() # at least 1 X
        if (self.Y_chromosome_template is None) or (sex == Genome.FEMALE):
            r = 0
        elif sex == Genome.MALE:
            r = 1
        else:
            r = random.randint(0,1)
        chrom_2_template = self.X_chromosome_template if r == 0 else self.Y_chromosome_template
        chrom_2 = chrom_2_template.generate_hom_recessive_chromosome()
        return ChromosomeSexPair([chrom_1,chrom_2], self.X_chromosome_template, chrom_2_template)

    def generate_het_sex_pair(self, sex=None, ploidy=None, rand_phase=None):

        if sex == Genome.MALE:
            # Can't generate heterozygous. Generate random alleles
            chrom_1_template = self.X_chromosome_template
            chrom_2_template = self.Y_chromosome_template
            chrom_1 = chrom_1_template.generate_random_chromosome()
            chrom_2 = chrom_2_template.generate_random_chromosome()
            return ChromosomeSexPair([chrom_1, chrom_2], self.X_chromosome_template, self.Y_chromosome_template)

        if self.X_chromosome_template is None:
            return None

        pair = self.X_chromosome_template.generate_het_pair(self.ploidy, rand_phase=rand_phase)
        return ChromosomeSexPair(pair.chrom_pair, self.X_chromosome_template, self.X_chromosome_template)


    def generate_random_sex_pair(self, sex=None):
        if self.X_chromosome_template is None:
            return None


        chrom_1 = self.X_chromosome_template.generate_random_chromosome() # at least 1 X
        if (self.Y_chromosome_template is None) or (sex == Genome.FEMALE):
            r = 0
        elif sex == Genome.MALE:
            r = 1
        else:
            r = random.randint(0,1)
        chrom_2_template = self.X_chromosome_template if r == 0 else self.Y_chromosome_template
        chrom_2 = chrom_2_template.generate_random_chromosome()
        return ChromosomeSexPair([chrom_1,chrom_2], self.X_chromosome_template, chrom_2_template)

    def generate_random_genome(self, sex=None):
        chrom_pairs = []
        for chromosome in self.chromosome_templates:
            chrom_pairs.append(chromosome.generate_random_pair(self.ploidy))
        sex_pair = self.generate_random_sex_pair(sex=sex)
        return Genome(self,chrom_pairs, sex_pair = sex_pair)

    def generate_hom_recessive_genome(self, sex=None):
        chrom_pairs = []
        for chromosome in self.chromosome_templates:
            chrom_pairs.append(chromosome.generate_hom_recessive_pair(self.ploidy))
        sex_pair = self.generate_hom_recessive_sex_pair(sex=sex)
        return Genome(self,chrom_pairs,sex_pair = sex_pair)

    def generate_hom_dominant_genome(self, sex=None):
        chrom_pairs = []
        for chromosome in self.chromosome_templates:
            chrom_pairs.append(chromosome.generate_hom_dominant_pair(self.ploidy))
        sex_pair = self.generate_hom_dominant_sex_pair(sex=sex)
        return Genome(self, chrom_pairs, sex_pair=sex_pair)

    def generate_het_genome(self, rand_phase=False, sex=None):
        chrom_pairs = []
        for chromosome in self.chromosome_templates:
            chrom_pairs.append(chromosome.generate_het_pair(self.ploidy,rand_phase=rand_phase))
        sex_pair = self.generate_het_sex_pair(sex=sex, rand_phase=rand_phase)
        return Genome(self,chrom_pairs, sex_pair=sex_pair)



class Chromosome(SerialiserMixin):
    def __init__(self,chromosome_template,alleles):
        self.chromosome_template = chromosome_template
        self.alleles = alleles

    # @staticmethod
    # def generate_random_chromosome(genes):
    #     actual_alleles = []
    #     for gene in genes:
    #         r = random.randint(0, len(gene) - 1)
    #         actual_alleles.append(gene.alleles[r])
    #     return Chromosome(actual_alleles)

    def __copy__(self):
        new_alleles = self.alleles.copy()
        return Chromosome(self.chromosome_template,new_alleles)

    def __str__(self):

       chrom_gen_str = ''.join([str(allele) for allele in self.alleles])
       return chrom_gen_str

    def alleles_order_by_lowest_alpha(self):
        if str(self.alleles[0]).upper() > str(self.alleles[-1]).upper():
            return ''.join([str(allele) for allele in self.alleles[::-1]])
        else:
            return ''.join([str(allele) for allele in self.alleles])



    @staticmethod
    def _from_attr_dict(attr_dict):
        alleles = [Allele._from_attr_dict(al) for al in attr_dict['alleles']]
        chrom_template = ChromosomeTemplate._from_attr_dict(attr_dict['chromosome_template'])
        obj = Chromosome(chrom_template,alleles)
        return obj



class ChromosomePair(SerialiserMixin):

    UNLINKED_THRESHOLD = 50000000
    CROSSOVER_DIVISOR  = 1000000

    def __init__(self,chromosome_template,chrom_pair):
        self.chromosome_template = chromosome_template
        self.chrom_pair = chrom_pair

    def __str__(self):
        out_str = ''
        for i in range(self.num_genes()):
            allele_pair = self.get_allele_pair(i)
            out_str +=  ''.join(allele_pair)

        return out_str

    @staticmethod
    def _from_attr_dict(attr_dict):
        chrom_pair =  [Chromosome._from_attr_dict(chrom) for chrom in attr_dict['chrom_pair']]
        chrom_template = ChromosomeTemplate._from_attr_dict(attr_dict['chromosome_template'])
        obj = ChromosomePair(chrom_template, chrom_pair)
        return obj

    def get_phase(self,alpha_sort=False, show_dist = True):
        alleles = []
        if alpha_sort:
            alleles = [chrom.alleles_order_by_lowest_alpha() for chrom in self.chrom_pair]
        else:
            alleles = [str(chrom) for chrom in self.chrom_pair]

        if alpha_sort:
           alleles.sort()

        if show_dist:
           prev_gen = None
           dists = []
           for i, allele_symbol in enumerate(alleles[0]):
               gene = self.chromosome_template.get_gene_from_symbol(allele_symbol)
               if i == 0:
                   prev_gene = gene
               else:
                   dists.append(gene.distance(prev_gene))
                   prev_gene = gene
           dists = [math.floor(abs(dist / 10000000)) * '-' for dist in dists]
           alleles_with_dist = ['','']
           for i in range(len(alleles[0])):
               if i == 0:
                  alleles_with_dist[0] = alleles[0][i]
                  alleles_with_dist[1] = alleles[1][i]
               else:
                  alleles_with_dist[0] += dists[i-1] + alleles[0][i]
                  alleles_with_dist[1] += dists[i - 1] + alleles[1][i]

           alleles = alleles_with_dist
           #alleles[0] = alleles[0][0] + dists[0] + alleles[0][1] + dists[1] + alleles[0][2]
           #alleles[1] = alleles[1][0] + dists[0] + alleles[1][1] + dists[1] + alleles[1][2]



        out_str = '//'.join(alleles)
        return out_str

    def get_parental_gamete(self,parent_num):
        gamete = ''
        for allele in self.chrom_pair[parent_num].alleles:
            gamete += allele.symbol
        return gamete




    def num_genes(self):
        return len(self.chromosome_template.genes)


    def get_allele_pair(self,i, sort=True):
        allele_pair= []
        for chrom in self.chrom_pair:
            allele_pair.append(str(chrom.alleles[i]))
        if sort:
            allele_pair.sort()
        return allele_pair

    def phenotype(self,alpha_sort=True):
        phen = ''
        allele_pairs = []
        for i in range(self.num_genes()):
            allele_pair = self.get_allele_pair(i)
            allele_pairs.append(allele_pair)

        if alpha_sort:
            allele_pairs.sort(key=lambda x:x[1])
        for allele_pair in allele_pairs:
            num_lower = 0
            lower_allele = ''

            for allele in allele_pair:
                if not (allele.isupper()):
                    num_lower +=1
                lower_allele = allele.lower()


            phen += lower_allele
            phen += '+' if num_lower < 2 else '-'
        if len(phen) == 6:
            pass
        elif (len(phen) == 4) :
            #phen += 'c+'
            pass
        elif (len(phen) == 2):
            #phen += 'b+c+'
            pass

        return phen

    def phenotype_afflicted(self,alpha_sort=True):
        afflicted_list = []
        for i, gene in enumerate(self.chromosome_template.genes):
            allele_pair = self.get_allele_pair(i)
            num_lower = 0
            for allele in allele_pair:
                if not (allele.isupper()):
                    num_lower += 1
            aff = ''
            if gene.inheritance_pattern == Gene.INH_PATTERN_DOMINANT:
               if num_lower > 0:
                  aff = allele_pair[0].lower()
               else:
                  aff =  allele_pair[0].upper()
            else:
                if num_lower ==  2:
                    aff = allele_pair[0].lower()
                else:
                    aff = allele_pair[0].upper()

            afflicted_list.append(aff)

        if alpha_sort:
            afflicted_list.sort(key=lambda x:x.upper())

        return ''.join(afflicted_list)

    def pick_random_alleles(self):

        new_alleles = []
        for i in range(self.num_genes()):
            r = random.randint(0, len(self.chrom_pair)-1)
            new_alleles.append(self.chrom_pair[r].alleles[i])

        return Chromosome(self.chromosome_template,new_alleles)


    def all_unlinked(self):
        all_positions = [gene.position for gene in self.chromosome_template.genes]
        highest = max(all_positions)
        lowest = min(all_positions)
        if ((highest - lowest) >=  ChromosomePair.UNLINKED_THRESHOLD):
            return True
        return False

    def allele_pairs(self):
        pairs = []
        for  i in range(self.num_genes()):
            pair = [(self.chrom_pair[0].alleles[i],0),(self.chrom_pair[1].alleles[i],1)]
            pairs.append(pair)

        return pairs

    def possible_gametes(self, suppress_combine_same=False):
        probs = self.crossover_probabilities()


        possibles = list(itertools.product(*self.allele_pairs()))

        possibles_with_prob = []

        for possible in possibles:
            prob = 0.5
            for i, (allele,phase) in enumerate(possible):
                if i == 0:
                    prev_phase = phase
                else:
                    this_prob = 1 - probs[i-1] if phase == prev_phase else  probs[i-1]
                    prob *= this_prob
                    prev_phase = phase


                #print(allele,phase)

            alleles_only = [p[0] for p in possible]
            if suppress_combine_same:
                prob = 1 / len(possibles)
                possibles_with_prob.append([alleles_only, prob, 1, len(possibles) ])
            else:
                found = False
                for p_with_prob in possibles_with_prob:
                    if alleles_only == p_with_prob[0]:
                        p_with_prob[1] += prob
                        p_with_prob[2] += 1
                        found = True
                        break
                if not found:
                    possibles_with_prob.append([alleles_only, prob, 1, len(possibles)])

        return possibles_with_prob

    def crossover_probabilities(self):

        probs = []

        for i in range(self.num_genes()-1):
            gene1 = self.chromosome_template.genes[i]
            gene2 = self.chromosome_template.genes[i+1]
            dist = gene1.distance(gene2)
            if debug > 0:
                print('dist: ', gene1, gene2, dist)
            prob_crossover = (dist / ChromosomePair.CROSSOVER_DIVISOR) / 100.0
            if debug > 0:
                print('prob: ',prob_crossover)
            if prob_crossover > 0.5:
                prob_crossover = 0.5
            probs.append(prob_crossover)

        return probs

    def meiosis(self):

        probs = self.crossover_probabilities()

        crossovers = []
        for i in range(self.num_genes()-1):
            gene1 = self.chromosome_template.genes[i]
            gene2 = self.chromosome_template.genes[i+1]
            dist = gene1.distance(gene2)
            if debug > 0:
                print('dist: ', gene1, gene2, dist)
            prob_crossover = (dist / ChromosomePair.CROSSOVER_DIVISOR) / 100.0
            if debug > 0:
                print('prob: ',prob_crossover)
            if prob_crossover >= 0.5:
                crossovers.append('Rand')
            else:
                r = random.random()
                if r < prob_crossover:
                    crossovers.append('Yes')
                else:
                    crossovers.append('No')

        if debug > 0:
            print(crossovers)


        if len(crossovers) == 0:  # Only 1 gene. Pick chromosome at random
            r = random.randint(0, len(self.chrom_pair) - 1)
            return self.chrom_pair[r]

        new_alleles = []
        for i,cr in enumerate(crossovers):
            if i == 0: # Pick which chrom pair to start with
                pair_num = random.randint(0, len(self.chrom_pair) - 1)
                if debug > 0:
                    print('first choice: ', pair_num, self.chrom_pair[pair_num].alleles[i])
                new_alleles.append(self.chrom_pair[pair_num].alleles[i])
            if crossovers[i] == 'Rand':
                pair_num = random.randint(0, len(self.chrom_pair) - 1)
                if debug > 0:
                    print('random choice', pair_num, self.chrom_pair[pair_num].alleles[i+1])
                new_alleles.append(self.chrom_pair[pair_num].alleles[i+1])
            else:
                if crossovers[i] == 'Yes':
                    pair_num = 0 if pair_num == 1 else 1
                if debug > 0:
                    print(crossovers[i] + ' choice',pair_num,self.chrom_pair[pair_num].alleles[i+1] )
                new_alleles.append(self.chrom_pair[pair_num].alleles[i+1])

        if debug > 0:
            print('new alleles: ')
            for allele in new_alleles:
                print(allele)

        return Chromosome(self.chromosome_template, new_alleles)



    def mate(self,other_chrom_pair):

        # if (self.num_genes() == 1):
        #     new_chrom_pair = []
        #     r = random.randint(0,len(self.chrom_pair)-1)
        #     new_chrom_pair.append(self.chrom_pair[r])
        #     r = random.randint(0,len(other_chrom_pair.chrom_pair)-1)
        #     new_chrom_pair.append(other_chrom_pair.chrom_pair[r])
        #
        #     return ChromosomePair(self.chromosome_template,new_chrom_pair)
        # else:
        #     return ChromosomePair(self.chromosome_template,[self.meiosis(),other_chrom_pair.meiosis()])

        return ChromosomePair(self.chromosome_template, [self.meiosis(), other_chrom_pair.meiosis()])

class ChromosomeSexPair(ChromosomePair, SerialiserMixin):
    def __init__(self,chrom_pair, chrom_1_template, chrom_2_template):
        if chrom_1_template.type == ChromosomeTemplate.X:
            super(ChromosomeSexPair,self).__init__(chrom_1_template, chrom_pair)
        else:
            super(ChromosomeSexPair, self).__init__(chrom_2_template, chrom_pair)

        self.chrom_1_template = chrom_1_template # ie first of chrom pair
        self.chrom_2_template = chrom_2_template # ie second of chrom pair

    def get_allele_X_pair(self,i, sort=True):
        pass

    def get_allele_Y(self,i, sort=True):
        pass

    def get_allele_pair(self, i, sort=True):
        # TODO fill in allele pair logic:
        allele_pair = []
        if self.chrom_1_template == self.chrom_2_template: #Assume female
           for chrom in self.chrom_pair:
               allele_pair.append('X' + str(chrom.alleles[i]))

        else:
           if self.chrom_1_template.type == ChromosomeTemplate.X:
               allele_pair.append('X' + str(self.chrom_pair[0].alleles[i]))
               allele_pair.append('Y_')

        if sort:
            allele_pair.sort()
        return allele_pair

    def get_allele_X_pair(self, i, sort=True):
            # TODO fill in allele pair logic:
            allele_pair = []
            if self.chrom_1_template == self.chrom_2_template:  # Assume female
                for chrom in self.chrom_pair:
                    allele_pair.append('X' + str(chrom.alleles[i]))

            else:
                if self.chrom_1_template.type == ChromosomeTemplate.X:
                    allele_pair.append('X' + str(self.chrom_pair[0].alleles[i]))
                    allele_pair.append('Y ')
                else:
                    allele_pair.append('X' + str(self.chrom_pair[1].alleles[i]))
                    allele_pair.append('Y ')

            if sort:
                allele_pair.sort()
            return allele_pair

    def get_allele_Y_pair(self, i, sort=True):
            allele_pair = []
            if self.chrom_1_template == self.chrom_2_template:
                return ''

            if self.chrom_1_template.type == ChromosomeTemplate.Y:
                allele_pair.append('Y' + str(self.chrom_pair[0].alleles[i]))
                allele_pair.append('X ')
            else:
                allele_pair.append('Y' + str(self.chrom_pair[1].alleles[i]))
                allele_pair.append('X ')


            if sort:
                allele_pair.sort()
            return allele_pair

        #self.get_allele_X_pair(i, sort=sort)
        # default to just X-linked

        # allele_pair = []
        # for chrom in self.chrom_pair:
        #     allele_pair.append(str(chrom.alleles[i]))
        # if sort:
        #     allele_pair.sort()
        # return allele_pair

    def phenotype_afflicted(self,alpha_sort=True):
        afflicted_list = []
        if self.chrom_1_template == self.chrom_2_template:
           x_chrom_template = self.chrom_1_template
           y_chrom_template = None
        else:
           if self.chrom_1_template.type == ChromosomeTemplate.Y:
               y_chrom_template = self.chrom_1_template
               x_chrom_template = self.chrom_2_template
           else:
               y_chrom_template = self.chrom_2_template
               x_chrom_template = self.chrom_1_template

        for i, gene in enumerate(x_chrom_template.genes):
            allele_pair = self.get_allele_X_pair(i)
            num_lower = 0
            num_upper = 0
            aff = ''
            for allele in allele_pair:
                if allele[0] == 'X': # female allele
                    if allele[1].isupper():
                        num_upper += 1
                    else:
                        num_lower +=1
            if y_chrom_template is not None:
                aff = allele_pair[0][1]
                # if allele_pair[0][1].isupper():
                #     aff = allele_pair[0][1].lower()
                # else:
                #     aff = allele_pair[0][1].upper()
            else:

                if gene.inheritance_pattern == Gene.INH_PATTERN_DOMINANT:
                   if num_lower > 0:
                      aff = allele_pair[0][1].lower()
                   else:
                      aff =  allele_pair[0][1].upper()
                else:
                    if num_upper > 0:
                        aff = allele_pair[0][1].upper()
                    else:
                        aff = allele_pair[0][1].lower()

            afflicted_list.append(aff)

        if y_chrom_template is None:
                #TODO Need to set afflicted to no for all  Y genes (even though this is a female and doesn't have a Y template)
                pass
        else:
            for i, gene in enumerate(y_chrom_template.genes):
                    allele_pair = self.get_allele_Y_pair(i)
                    num_lower = 0
                    num_upper = 0
                    aff = ''
                    for allele in allele_pair:
                        if allele[0] == 'X': # female allele
                            pass
                        else:
                            if allele[1].isupper():
                                num_upper += 1
                                aff = allele[1]
                            else:
                                num_lower +=1
                                aff = allele[1]


                    # if gene.inheritance_pattern == Gene.INH_PATTERN_DOMINANT:
                    #    if num_upper >= 1:
                    #       aff = allele_pair[0][1].upper()
                    #    else:
                    #       aff =  allele_pair[0][1].lower()
                    # else:
                    #     if num_upper >= 1:
                    #         aff = allele_pair[0][1].lower()
                    #     else:
                    #         aff = allele_pair[0][1].upper()

                    afflicted_list.append(aff)



        if alpha_sort:
            afflicted_list.sort(key=lambda x:x.upper())

        return ''.join(afflicted_list)


    def num_genes(self):
         return  self.num_X_genes() + self.num_Y_genes() #len(self.chromosome_template.genes)

    def num_X_genes(self):
        if self.chrom_1_template.type == ChromosomeTemplate.X:
            return len(self.chrom_1_template.genes)
        else:
            return len(self.chrom_2_template.genes)

    def num_Y_genes(self):
        if self.chrom_1_template == self.chrom_2_template:
            return 0
        if self.chrom_1_template.type == ChromosomeTemplate.Y:
            return len(self.chrom_1_template.genes)
        else:
            return len(self.chrom_2_template.genes)

    def meiosis(self):

        if self.chrom_1_template == self.chrom_2_template: #assume female
            return super().meiosis()
        else:
            r = random.randint(0,1) # male or female!
            return self.chrom_pair[r]

        # crossovers = []
        # for i in range(self.num_genes() - 1):
        #     gene1 = self.chromosome_template.genes[i]
        #     gene2 = self.chromosome_template.genes[i + 1]
        #     dist = gene1.distance(gene2)
        #     if debug > 0:
        #         print('dist: ', gene1, gene2, dist)
        #     prob_crossover = (dist / 1000000) / 100.0
        #     if debug > 0:
        #         print('prob: ', prob_crossover)
        #     if prob_crossover >= 0.5:
        #         crossovers.append('Rand')
        #     else:
        #         r = random.random()
        #         if r < prob_crossover:
        #             crossovers.append('Yes')
        #         else:
        #             crossovers.append('No')
        #
        # if debug > 0:
        #     print(crossovers)
        #
        # if len(crossovers) == 0:  # Only 1 gene. Pick chromosome at random
        #     r = random.randint(0, len(self.chrom_pair) - 1)
        #     return self.chrom_pair[r]
        #
        # new_alleles = []
        # for i, cr in enumerate(crossovers):
        #     if i == 0:  # Pick which chrom pair to start with
        #         pair_num = random.randint(0, len(self.chrom_pair) - 1)
        #         if debug > 0:
        #             print('first choice: ', pair_num, self.chrom_pair[pair_num].alleles[i])
        #         new_alleles.append(self.chrom_pair[pair_num].alleles[i])
        #     if crossovers[i] == 'Rand':
        #         pair_num = random.randint(0, len(self.chrom_pair) - 1)
        #         if debug > 0:
        #             print('random choice', pair_num, self.chrom_pair[pair_num].alleles[i + 1])
        #         new_alleles.append(self.chrom_pair[pair_num].alleles[i + 1])
        #     else:
        #         if crossovers[i] == 'Yes':
        #             pair_num = 0 if pair_num == 1 else 1
        #         if debug > 0:
        #             print(crossovers[i] + ' choice', pair_num, self.chrom_pair[pair_num].alleles[i + 1])
        #         new_alleles.append(self.chrom_pair[pair_num].alleles[i + 1])
        #
        # if debug > 0:
        #     print('new alleles: ')
        #     for allele in new_alleles:
        #         print(allele)
        #
        # return Chromosome(self.chromosome_template, new_alleles)

    def mate(self, other_chrom_pair):

        # if (self.num_genes() == 1):
        #     new_chrom_pair = []
        #     r = random.randint(0,len(self.chrom_pair)-1)
        #     new_chrom_pair.append(self.chrom_pair[r])
        #     r = random.randint(0,len(other_chrom_pair.chrom_pair)-1)
        #     new_chrom_pair.append(other_chrom_pair.chrom_pair[r])
        #
        #     return ChromosomePair(self.chromosome_template,new_chrom_pair)
        # else:
        #     return ChromosomePair(self.chromosome_template,[self.meiosis(),other_chrom_pair.meiosis()])

        chrom_1 = self.meiosis()
        chrom_2 = other_chrom_pair.meiosis()

        return ChromosomeSexPair([chrom_1,chrom_2], chrom_1.chromosome_template, chrom_2.chromosome_template)




    def __str__(self):
        out_str = ''

        for i in range(self.num_X_genes()):
            allele_pair = self.get_allele_X_pair(i)
            out_str +=  ''.join(allele_pair)

        out_str += ' '

        for i in range(self.num_Y_genes()):
            allele_pair = self.get_allele_Y_pair(i)
            out_str += ''.join(allele_pair)

        return out_str


class Genome(SerialiserMixin):

    MALE = 'M'
    FEMALE = 'F'
    UNKNOWN = 'U'

    def __init__(self,genome_template,chromosome_pairs=[],sex_pair=None):
        self.genome_template = genome_template
        self.chromosome_pairs = chromosome_pairs
        self.sex_pair =  sex_pair

    def __str__(self):

        out_str = ''

        #out_str += self.sex() + ': '
        for chrom_pair in self.chromosome_pairs:
            out_str += str(chrom_pair)



        str_list = []
        for i in range(int(len(out_str) / 2)):
            str_list.append(out_str[i * 2:i * 2 + 2])
        str_list.sort(key=lambda x: x[0].upper())

        out_str = ''.join(str_list)

        if self.sex_pair is None:
            pass
        else:
            out_str += '  ' + str(self.sex_pair)

        if self.sex == Genome.UNKNOWN:
            return out_str
        else:
            return self.sex + ': ' + out_str
        #return out_str

    @staticmethod
    def _from_attr_dict(attr_dict):
        chrom_pairs = [ChromosomePair._from_attr_dict(chrom_pair) for chrom_pair in attr_dict['chromosome_pairs']]
        genome_template =GenomeTemplate._from_attr_dict(attr_dict['genome_template'])
        obj = Genome(genome_template,chromosome_pairs=chrom_pairs)
        return obj

    def get_phase(self,alpha_sort=False):
        chrom_phases = [chrom_pair.get_phase(alpha_sort) for chrom_pair in self.chromosome_pairs]
        if alpha_sort:
           chrom_phases.sort(key=lambda x: x.upper())
        phase = ';'.join(chrom_phases)
        return phase

    def get_parental_gamete(self,parent_num,sort_alpha=False):
        gamete_per_pair =  [chromosome_pair.get_parental_gamete(parent_num) for chromosome_pair in self.chromosome_pairs]
        gamete = ''.join(gamete_per_pair)
        if sort_alpha:
           gamete_chr_list = list(gamete)
           gamete_chr_list.sort(key=lambda x: x.upper())
           gamete = ''.join(gamete_chr_list)
        return gamete

    def genotype(self):
        p1 = []
        p2 = []
        for  chrom_pair in self.chromosome_pairs:
            p1.append(str(chrom_pair.chrom_pair[0]))
            p2.append(str(chrom_pair.chrom_pair[1]))
        if self.sex_pair is None:
            pass
        else:
            p1.append(str(self.sex_pair.chrom_pair[0]))
            p2.append(str(self.sex_pair.chrom_pair[1]))

        return '\n'.join(['//'.join(p1),'//'.join(p2)])

    def phenotype(self,sort_alpha=True):
        phen = ''
        for chromosome_pair in self.chromosome_pairs:
            phen += chromosome_pair.phenotype()

        if not 'a' in phen:
           phen += 'a+'
        if not 'b' in phen:
           phen += 'b+'
        if not 'c' in phen:
           phen += 'c+'

        if sort_alpha:
            phen_list = []
            for i in range(int(len(phen) / 2)):
                phen_list.append(phen[i*2:i*2+2])
            phen_list.sort(key=lambda x: x[0])
            return ''.join(phen_list)


        return phen

    def phenotype_afflicted(self,sort_alpha=True):
        aff = ''
        for chromosome_pair in self.chromosome_pairs:
            aff += chromosome_pair.phenotype_afflicted()

        if self.sex_pair is None:
            pass
        else:
            aff += self.sex_pair.phenotype_afflicted()

        if sort_alpha:
            aff_list = list(aff)
            aff_list.sort(key=lambda x: x.upper())
            aff = ''.join(aff_list)

        return aff


    def gen_phen(self):
        return {'gen':str(self),'phen':self.phenotype(), 'gen_phase':self.genotype()}


    def possible_gametes(self, suppress_combine_same=False):
        all_possibles_with_prob = []
        for i, cp in enumerate(self.chromosome_pairs):
            possibles_with_prob  = cp.possible_gametes(suppress_combine_same=suppress_combine_same)
            all_possibles_with_prob.append(possibles_with_prob)

        possibles = list(itertools.product(*all_possibles_with_prob))

        combined_possibles = []
        for possible in possibles:
            combined_alleles = []
            combined_freq = 1
            combined_prob = 1
            combined_tot_possible = 1
            for chrom_possible in possible:
                combined_alleles.extend(chrom_possible[0])
                combined_freq *= chrom_possible[1]
                combined_prob *= chrom_possible[2]
                combined_tot_possible *= chrom_possible[3]
            combined_possibles.append([combined_alleles, combined_freq, combined_prob, combined_tot_possible])


        return combined_possibles

    def possible_gametes_formatted(self, dec_places=3, suppress_combine_same = False):
        possibles = [[''.join([str(p) for p in possible[0]]),round(possible[1],dec_places),possible[2],possible[3]] for possible in self.possible_gametes(suppress_combine_same=suppress_combine_same)]

        possibles.sort(key=lambda x: x[0])
        return possibles

    def mate(self, other_genome):
        new_chrom_pairs = []
        for i,cp in enumerate(self.chromosome_pairs):
            new_chrom_pairs.append(cp.mate(other_genome.chromosome_pairs[i]))

        if self.sex_pair is None:
            new_sex_pair = None
        else:
            new_sex_pair = self.sex_pair.mate(other_genome.sex_pair)
        return Genome(self.genome_template,new_chrom_pairs, new_sex_pair)

    @property
    def sex(self):
        if self.sex_pair is None:
            return Genome.UNKNOWN

        for chrom in self.sex_pair.chrom_pair:
            if chrom.chromosome_template.type == ChromosomeTemplate.Y:
                return Genome.MALE
        return Genome.FEMALE

    @staticmethod
    def test_cross_het_gametes_to_phenotypes(het_gametes=['ABC','ABc', 'AbC', 'Abc', 'aBC', 'aBc', 'abC', 'abc']):
       phenotypes = {gamete: ''.join([ch + '+' if ch.isupper() else ch + '-' for ch in gamete]) for gamete in het_gametes}
       return phenotypes


class Organism(SerialiserMixin):

    def __init__(self,genome, partner=None, parents=None, children=None, counter_id = None, level=None):
       self.genome = genome
       if parents is None:
           self._parents  = []
       else:
           self._parents = parents
       if level is None:
           if self.has_parents:
               self.level = self.parents[0].level + 1
           else:
               self.level = 0
       else:
           self.level = level
       self.set_partner(partner)
       if children is None:
          self._children = []
       else:
          self._children = children
       if counter_id is None:
          self.counter_id = 1
       else:
          self.counter_id = counter_id

       self._possible_genotypes = {} #Inferrable genotypes for each consistent pedigree


    @staticmethod
    def organism_with_random_genotype(genome_template, sex=None, counter_id=None, level=None):
        genome = genome_template.generate_random_genome(sex=sex)
        return Organism(genome, counter_id=counter_id, level=level)

    @staticmethod
    def organism_with_hom_recessive_genotype(genome_template, sex=None, counter_id=None, level=None):
        genome = genome_template.generate_hom_recessive_genome(sex=sex)
        return Organism(genome, counter_id=counter_id, level=level)

    @staticmethod
    def organism_with_hom_dominant_genotype(genome_template, sex=None, counter_id=None, level=None):
        genome = genome_template.generate_hom_dominant_genome(sex=sex)
        return Organism(genome, counter_id=counter_id, level=level)

    @staticmethod
    def organism_with_het_genotype(genome_template,rand_phase=False, sex=None, counter_id=None, level=None):
        genome = genome_template.generate_het_genome(rand_phase=rand_phase, sex=sex)
        return Organism(genome, counter_id=counter_id, level=level)

    # @staticmethod
    # def generate_pedigree(max_levels=2, type='auto_rec', prob_mate = 0.5, max_children=4):
    #     if type != 'auto_rec':
    #         raise Exception('pedigree type not allowed: ' + type)
    #
    #     g = Gene(AlleleSet.default_alleleset_from_symbol('A'), 1000000,
    #                   inheritance_pattern=Gene.INH_PATTERN_RECESSIVE)
    #
    #     c = ChromosomeTemplate('PedAutoChrom', 2000000, [g1])
    #     c_X = ChromosomeTemplate('PedXChrom', 2000000, [], type=ChromosomeTemplate.X)
    #     c_Y = ChromosomeTemplate('PedYChrom', 2000000, [], type=ChromosomeTemplate.Y)
    #
    #     gt = GenomeTemplate(ploidy=2, chromosome_templates=[c], X_chromosome_template=c_X,
    #                         Y_chromosome_template=c_Y, name='PedGT')
    #     print(str(gt))
    #
    #     next_id = 1
    #     adam = Organism.organism_with_random_genotype(gt, sex=Genome.MALE, counter_id=next_id)
    #     print(str(adam))
    #     next_id += 1
    #     # eve = Organism.organism_with_random_genotype(gt, sex=Genome.FEMALE, counter_id=next_id)
    #     # print(str(eve))
    #     # next_id += 1
    #
    #     for i in range(0,max_levels-1):
    #         orgs = adam.all_orgs_in_pedigree(include_in_laws=False, level=i)
    #         for org in orgs:
    #             r = random.randint(0,1)
    #             if i == 0 or r <=  prob_mate:
    #                 org_mate = Organism.organism_with_random_genotype(gt,
    #                                                                   sex=Genome.MALE if org.sex == Genome.FEMALE else Genome.FEMALE,
    #                                                                   counter_id=next_id)
    #                 next_id += 1
    #                 org.set_partner(org_mate)
    #                 num_ch = random.randint(0,max_children)
    #                 org.mate(times=num_ch, next_id=next_id)
    #                 next_id += org.num_children
    #
    #     return adam

    @staticmethod
    def calc_recombination_fractions(gametes):
        sorted_gametes = sorted(gametes.items(), key=lambda x: x[1], reverse=True)
        print(sorted_gametes)
        parentals = [sorted_gametes[0][0],sorted_gametes[1][0] ]
        double_recombinations = [sorted_gametes[-2][0],sorted_gametes[-1][0] ]
        print('parentals', parentals)
        print('double recombs',double_recombinations)

        diffs = []
        sames = []
        for i, ch in enumerate(parentals[0]):
            if ch == double_recombinations[0][i]:
               sames.append(ch)
            else:
               diffs.append(ch)
        middle_gene = sames[0] if len(sames) == 1 else diffs[0]
        print('middle: ', middle_gene)

        outer_genes = []
        for i, ch in enumerate(parentals[0]):
             if ch == middle_gene:
                 pass
             else:
                 outer_genes.append(ch)
        parental_ordered = outer_genes[0] + middle_gene + outer_genes[1]
        print('parental ordered: ', parental_ordered)
        parental_other_ordered = parental_ordered.swapcase()
        print('parental other ordered: ', parental_other_ordered)

        def get_pair(item, ind1, ind2):
            return item[ind1] + item[ind2]



        # for phen_combs in phen_combinations_per_pair:
        #     observed = [count for key, count in phen_combs.items()]
        #     chisq, p = chisquare(observed, ddof=1)
        #     phen_combs['abbrev'] = list(phen_combs.keys())[0].replace('+', '')
        #     phen_combs['abbrev'] = phen_combs['abbrev'].replace('-', '')
        #     phen_combs['p'] = '{:.2e}'.format(p)

        recombination_counts = [0,0,0]
        tot_gametes = 0

        for gamete, count in gametes.items():
            tot_gametes+= count
            all_pair_inds = [[0,1],[1,2],[0,2]]
            for i, pair_inds in enumerate(all_pair_inds):
                gamete_pair = get_pair(gamete, pair_inds[0],pair_inds[1])
                parental_pair = get_pair(parentals[0], pair_inds[0],pair_inds[1])
                if (gamete_pair == parental_pair) or (gamete_pair == parental_pair.swapcase()):
                    pass
                else:
                    recombination_counts[i] += count

        p_values = []
        for count in recombination_counts:
            chisq, p = chisquare([count, tot_gametes - count], ddof = 0)
            p_values.append(p)

        recombination_pairs = ['ab', 'bc', 'ac']

        recombination_fractions = [[recombination_pairs[i],count / tot_gametes, "{0:.3e}".format(p_values[i])] for i, count in enumerate(recombination_counts)]


        recombination_indexes = [recombination_pairs.index(''.join(sorted(parental_ordered[:2].lower()))), recombination_pairs.index(''.join(sorted(parental_ordered[1:3].lower()) ))]
        linkages = ['L' if p_values[i] < 0.04 else 'U' for i in recombination_indexes]

        if linkages[0] == 'U' and linkages[1] == 'U':
            parental_ordered = 'ABC'
            parental_ordered_formatted = ';'.join([ch.upper() + '//' + ch.lower() for ch in parental_ordered])
        elif  linkages[0] == 'U' and linkages[1] == 'L':
            linked_part = parental_ordered[1:3]
            linked_part_list = list(linked_part)
            linked_part_list.sort(key=lambda x: x.upper())
            linked_part = ''.join(linked_part_list)
            if linked_part[0].islower():
                linked_part = linked_part.swapcase()
            unlinked_part = parental_ordered[0]
            if unlinked_part.islower():
                unlinked_part = unlinked_part.swapcase()
            parental_ordered = unlinked_part + linked_part
            parental_ordered_formatted = ';'.join([unlinked_part + '//' + unlinked_part.swapcase(),linked_part + '//' + linked_part.swapcase()])
        elif   linkages[0] == 'L' and linkages[1] == 'U':
            linked_part = parental_ordered[0:2]
            linked_part_list = list(linked_part)
            linked_part_list.sort(key=lambda x: x.upper())
            linked_part = ''.join(linked_part_list)
            if linked_part[0].islower():
                linked_part = linked_part.swapcase()
            unlinked_part = parental_ordered[2]
            if unlinked_part.islower():
                unlinked_part = unlinked_part.swapcase()
            parental_ordered = linked_part + unlinked_part
            parental_ordered_formatted = ';'.join([ linked_part + '//' + linked_part.swapcase(),unlinked_part + '//' + unlinked_part.swapcase()])
        else:
            if (parental_ordered[0].upper() > parental_ordered[2].upper()):
                parental_ordered = parental_ordered[2] + parental_ordered[1] + parental_ordered[0]
            if parental_ordered[0].islower():
               parental_ordered = parental_ordered.swapcase()
            parental_ordered_formatted = parental_ordered + '//' + parental_ordered.swapcase()

        for linkage in linkages:
            if linkage == 'U':
                pass


        return {'parentals': parentals, 'double_recombinations': double_recombinations, 'parental_ordered':parental_ordered, 'parental_ordered_formatted': parental_ordered_formatted, 'middle_gene': middle_gene,
                'recombination_fractions': recombination_fractions, 'p_values': p_values, 'linkages': ''.join(linkages), 'phenotypes': Genome.test_cross_het_gametes_to_phenotypes(gametes.keys())}

    @staticmethod
    def organism_from_gametes(gametes):
        fractions = Organism.calc_recombination_fractions(gametes)

        #genome = genome_template.generate_genome_from_gametes(gametes)
        #return Organism(genome)
        

    def to_json(self):
        """ Serialise to json for sending to javascript"""
        org = { 'id': self.counter_id,
                'sex': 'male' if self.sex == Genome.MALE else 'female',
                'afflicted': True if self.afflicted == 'a' else False,
                'level': self.level + 1, #javascript levels are 1 based
                'isInLaw': self.is_inlaw,
                'partner': None if self.partner is None else self.partner.counter_id,
                'children': [child.counter_id for child in self.children],
                'inferrable_genotypes': self._possible_genotypes
        }
        return org

    def mate(self,other_org=None, times=1, next_id=None):
        partner = None
        child = None

        if other_org is None:
            if self.partner is None:
                return None
            else:
                partner = self.partner

        else:
            partner = other_org

        for i in range(times):
            new_genome = self.genome.mate(partner.genome)

            if self.genome.sex == Genome.MALE:
               parent_1 = self
               parent_2 = partner
            else:
               parent_1  = partner
               parent_2 = self

            if next_id is None:
               counter_id = None
            else:
               counter_id = next_id + i
            child =  Organism(new_genome, parents = (parent_1, parent_2), counter_id=counter_id )
            parent_1.add_child(child)
            parent_2.add_child(child)

        return child

    def genotype(self):
        return self.genome.genotype()

    def gen_phen(self):
        return self.genome.gen_phen()

    def set_partner(self, partner):
        self.partner = partner
        if partner is None:
            pass
        else:
            partner.partner = self
            if partner.level > 0:
               self.level= partner.level
            elif self.level > 0:
                partner.level = self.level


    def set_parents(self, male_parent, female_parent):
        self._parents = (male_parent, female_parent)

    def add_child(self, child):
        self._children.append(child)

    def set_possible_genotype(self, inh_type, chrom_type, possible_gen):
        ''' Inferrable genotype'''
        key = chrom_type + str(inh_type)

        gen_set = False
        err = False

        if key in self._possible_genotypes:
           if self._possible_genotypes[key] == possible_gen:
               pass
           else:
               if self._possible_genotypes[key] is not None:
                   err = True
               self._possible_genotypes[key] = possible_gen
               gen_set = True

        else:
           self._possible_genotypes[key] = possible_gen
           gen_set = True

        return gen_set, err

    @property
    def is_inlaw(self):
        if self.level == 0:
            return False

        return not self.has_parents

    @property
    def children(self):
        return self._children

    @property
    def num_children(self):
        return len(self._children)

    @property
    def has_parents(self):
        return len(self._parents) > 0

    @property
    def parents(self):
        return self._parents

    @property
    def female_parent(self):
        for parent in self._parents:
            if parent.sex == Genome.FEMALE:
                return parent

        return None

    @property
    def male_parent(self):
        for parent in self._parents:
            if parent.sex == Genome.MALE:
                return parent

        return None

    @property
    def siblings(self):
        if not self.has_parents:
            return []

        _siblings = []
        for child in self.parents[0].children:
            if child is self:
                pass
            else:
                _siblings.append(child)
        return _siblings

    @property
    def next_younger_sibling(self):
        siblings = self.siblings
        lowest_younger_id = math.inf
        lowest_younger = None
        for sibling in siblings:
            if (sibling.counter_id > self.counter_id) and (sibling.counter_id < lowest_younger_id):
                lowest_younger = sibling
                lowest_younger_id = sibling.counter_id
        return lowest_younger

    @property
    def next_older_sibling(self):
        siblings = self.siblings
        highest_older_id = -1 * math.inf
        highest_older = None
        for sibling in siblings:
            if (sibling.counter_id < self.counter_id) and (sibling.counter_id > highest_older_id):
                highest_older = sibling
                highest_older_id = sibling.counter_id
        return highest_older

    @property
    def num_siblings(self):
        return len(self.siblings)

    def _get_nephews(self, include='all'):
        _nephews = []
        if include == 'all':
            for sibling in self.siblings:
                if sibling is self:
                    pass
                else:
                    _nephews.extend(sibling.children)

        elif include == 'o':
            next_older = self.next_older_sibling
            if next_older is None:
                pass
            else:
                _nephews = next_older.children
        elif include == 'y':
            next_younger = self.next_younger_sibling
            if next_younger is None:
                pass
            else:
                _nephews = next_younger.children

        return _nephews

    @property
    def next_older_nephews(self):
        return self._get_nephews(include='o')

    @property
    def next_younger_nephews(self):
        return self._get_nephews(include='y')

    @property
    def nephews(self):
        return self._get_nephews(include='all')

    @property
    def cousins(self):
        if self.parents is None:
            return []

        _cousins = []
        for parent in self.parents:
            for parent_sibling in parent.siblings:
                _cousins.extend(parent_sibling.children)
        return _cousins

    @property
    def num_cousins(self):
        return len(self.cousins)

    @property
    def sex(self):
        return self.genome.sex


    def orgs_afflicted_below(self, allele='a'):
        affl_list = []
        if allele in self.genome.phenotype_afflicted():
            affl_list.append(self.counter_id)
        if self.partner is None:
            pass
        else:
            if allele in self.partner.genome.phenotype_afflicted():
                affl_list.append(self.partner.counter_id)

        if self.num_children == 0:
            return affl_list
        else:
            affl_list_children = []
            for child in self.children:
                affl_list_children.extend(child.orgs_afflicted_below(allele=allele))
            affl_list.extend(affl_list_children)
            return affl_list


    # def peers_at_level(self,level,include_in_laws=False):
    #     pass

    # def all_orgs_in_pedigree(self, include_in_laws=False, level=None):
    #     adam = self.adam
    #     all_related = [adam]
    #     if adam.partner is not None:
    #         all_related.append(adam.partner)
    #     all_related.extend(adam.descendants)
    #     if include_in_laws:
    #         all_related.extend(adam.descendant_partners)
    #
    #     if level is None:
    #         return all_related
    #     else:
    #         all_related_with_level = [org for org in all_related if org.level == level]
    #         return all_related_with_level

    def print_org_tree_below(self):

        print('\t' * self.level + str(self) + ' affl: ' + str(self.genome.phenotype_afflicted()))
        if self.partner is None:
            pass
        else:
            print('\t' * self.partner.level  +'P '+ str(self.partner) + ' affl: ' + str(self.partner.genome.phenotype_afflicted()))
        for child in self.children:
            #print('\t' * child.level + str(child) + ' affl: ' + str(child.genome.phenotype_afflicted()))
            child.print_org_tree_below()



    @property
    def descendants(self):
        _descendants = []
        for child in self.children:
            _descendants.append(child)
            _descendants.extend(child.descendants)
        return _descendants

    @property
    def descendant_partners(self):
        _descendant_partners = []
        for child in self.children:
            if child.partner is None:
                pass
            else:
                _descendant_partners.append(child.partner)
            _descendant_partners.extend(child.descendant_partners)
        return _descendant_partners

    @property
    def num_descendant_partners(self):
        return len(self.descendant_partners)


    # @property
    # def adam(self):
    #     if self.parents is None:
    #        if self.sex == Genome.MALE:
    #             return self
    #        else:
    #            if self.partner is None:
    #                 return None
    #            else:
    #                 return self.partner.adam
    #
    #     org = self
    #     while org.parents is not None:
    #         org = org.parents[0]
    #     return org

    @property
    def num_descendants(self):
        return len(self.descendants)

    # @property
    # def total_organisms(self):
    #     if self.adam is None:
    #        return 1 if self.partner is None else 2
    #
    #     tot_orgs = self.adam.num_descendants + self.adam.num_descendant_partners +  1  # +1 for adam himself
    #     if self.adam.partner is None:
    #        pass
    #     else:
    #       tot_orgs +=1 # +1 for adams partner
    #     return tot_orgs

    @property
    def afflicted(self):
        return self.genome.phenotype_afflicted()

    def __str__(self):
        if self.counter_id is None:
            return str(self.genome)
        else:
            return f'{self.counter_id} - {self.genome}'

    @staticmethod
    def _from_attr_dict(attr_dict):
        ''' Note - serialisation currently only possible for Organism objects which have no parents or children (and hence have no recursion involved)'''
        genome =  Genome._from_attr_dict(attr_dict['genome'])
        obj = Organism(genome)
        return obj


    @staticmethod
    def unique_genotypes(organisms):
        unique_gens = {}
        for organism in organisms:
            if str(organism) in unique_gens:
                unique_gens[str(organism)] +=1
            else:
                unique_gens[str(organism)] = 1
        return unique_gens

    @staticmethod
    def unique_phenotypes(organisms):
        unique_phens = {}
        for organism in organisms:
            phen = organism.genome.phenotype()
            if phen in unique_phens:
                unique_phens[phen] +=1
            else:
                unique_phens[phen] = 1
        return unique_phens

class Pedigree:
    def __init__(self, max_levels=2, inh_type=Gene.INH_PATTERN_RECESSIVE, chrom_type = ChromosomeTemplate.AUTOSOMAL, symbol='A', prob_mate=0.5, max_children=4):
        self.max_levels = max_levels
        self.inh_type = inh_type
        self.chrom_type = chrom_type
        self.symbol = symbol
        self.prob_mate  = prob_mate
        self.max_children = max_children
        self.gt = None
        self.next_id = 1
        self._organisms = []
        self.gt = self.initialise_genome()

    @property
    def adam(self):

        for org in self._organisms:
            if org.parents is None and org.level == 0 and org.sex == Genome.MALE:
                return org

        return None

    @property
    def total_organisms(self):
        return len(self._organisms)

    def org_with_id(self, id):
        for org in self._organisms:
            if org.counter_id == id:
               return org
        return None

    def add_organism(self, sex, hom_rec=False, het=False, id=None, level=None):
        if id is None:
           id_to_allocate = self.next_id
           self.next_id +=1
        else:
           id_to_allocate = id
        if hom_rec:
            org = Organism.organism_with_hom_recessive_genotype(self.gt, sex=sex, counter_id=id_to_allocate, level=level)
        elif het:
            org = Organism.organism_with_het_genotype(self.gt, sex=sex, counter_id=id_to_allocate, level=level)
        else:
            org = Organism.organism_with_random_genotype(self.gt, sex=sex, counter_id=id_to_allocate, level=level)
        #print(str(org))
        self._organisms.append(org)
        #self.next_id += 1
        return org

    def all_orgs_in_pedigree(self, include_inlaws=True, level=None):

        if include_inlaws:
           orgs = self._organisms
        else:
           orgs= [org for org in self._organisms if not org.is_inlaw]

        if level is None:
            return orgs

        all_orgs_with_level = [org for org in orgs if org.level == level]
        return all_orgs_with_level

    def possible_genotypes(self, chrom_type=None, inh_type = None):

        if chrom_type is None:
            chrom_type =self.chrom_type
        if inh_type is None:
            inh_type = self.inh_type

        poss_gens = []
        done = False
        while not done:

            done = True

    def to_json(self):
        ped = { 'adam': 1,
                'orgs' : [org.to_json() for org in self._organisms]}

        return ped

    @staticmethod
    def pedigree_from_json(j):
       x = 1
       p = Pedigree(inh_type=Gene.INH_PATTERN_RECESSIVE, chrom_type=ChromosomeTemplate.AUTOSOMAL)

       def find_j_org_with_id_in_json(j, id):
           for j_org in j['orgs']:
               if j_org['id'] == id:
                   return j_org
           return None

       def add_org_for_j_org(j_org):
           if j_org['afflicted']:
               hom_rec = True
               het = False
           else:
               hom_rec = False
               het = True
           org = p.add_organism(sex='M' if j_org['sex'] == 'male' else 'F', hom_rec=hom_rec, het=het, id=j_org['id'], level=j_org['level'])
           return org

       for j_org in j['orgs']:
           org = p.org_with_id(j_org['id'])
           if org is not None:
               continue
           org = add_org_for_j_org(j_org) #p.add_organism(sex='M' if j_org['sex'] == 'male' else 'F', hom_rec = j_org['afflicted'], id=j_org['id'])
           org_mate = p.org_with_id(j_org['partner'])
           if org_mate is None:
              org_mate = add_org_for_j_org(find_j_org_with_id_in_json(j, j_org['partner']))
           org.set_partner(org_mate)

           for j_child_org_id in j_org['children']:

               child_org = add_org_for_j_org(find_j_org_with_id_in_json(j, j_child_org_id))
               male_parent = org if org.sex == Genome.MALE else org_mate
               female_parent = org if org.sex == Genome.FEMALE else org_mate
               child_org.set_parents(male_parent, female_parent)
               org.add_child(child_org)
               org_mate.add_child(child_org)

       act_gens = []
       for org in p.all_orgs_in_pedigree():
           act_gens.append(f'{org}')

       consistent_per_inferrer = {}
       possible_genotypes_per_inferrer = {}

       inferrers = [ARGenotypeInferrer(p), ADGenotypeInferrer(p), XRGenotypeInferrer(p), XDGenotypeInferrer(p),
                     YGenotypeInferrer(p)]
       for inferrer in inferrers:
            consistent, err_msg = inferrer.infer()

            if consistent:
                if inferrer.inferrer_type in consistent_per_inferrer:
                    consistent_per_inferrer[inferrer.inferrer_type] += 1
                else:
                    consistent_per_inferrer[inferrer.inferrer_type] = 1
            else:
                consistent_per_inferrer[inferrer.inferrer_type] = 0
            #     print(f'Consistent with {inferrer.inferrer_type}')

            possible_genotypes_per_inferrer[inferrer.inferrer_type] = inferrer.all_possible_genotypes


       return p, act_gens, consistent_per_inferrer, possible_genotypes_per_inferrer

    @staticmethod
    def pedigree_from_text(text):
        orgs = []
        ids_found = {}

        def convert_to_json_format(ids_found):
            j = {}
            for org_id, org in ids_found.items():
                if (org['level'] == 1) and (org['sex'] == 'male'):
                    j['adam'] = org_id
                    break
            orgs = []
            for org_id, org in ids_found.items():
                org_dict = {'id':org_id}
                for key, attr in org.items():
                    org_dict[key] = attr
                orgs.append(org_dict)
            j['orgs'] = orgs
            return j



        def register_id(org, ids_found, level=None, adding_child = False):
            m = re.match(r'(\d+)(.*)', org)
            id = int(m.group(1))
            attrs = m.group(2)
            if id not in ids_found:
               ids_found[id] = {'sex': 'male' if attrs[0] == 'M' else 'female', 'afflicted': True if (len(attrs) == 2 and attrs[1] == 'A') else False, 'partner' : None, 'inferrable_genotypes': [], 'children':[], 'isInLaw': False, 'level': level}
               if (not adding_child) and (level > 0):
                  ids_found[id]['isInLaw'] = True
            else:
                if level is not None:
                    ids_found[id]['level'] = level
            return id

        def parse_pedigree_text_line(line, ids_found, level=None):
            m = re.match(r'(.*)-(.*):(.*)', line)

            parent1 = m.group(1)
            p1_id = register_id(parent1, ids_found, level)
            if level is None:
                level = ids_found[p1_id]['level']
            parent2 = m.group(2)
            p2_id = register_id(parent2, ids_found, level)
            ids_found[p1_id]['partner'] = p2_id
            ids_found[p2_id]['partner'] = p1_id

            children = m.group(3).split(',')
            child_ids = []
            for child in children:
                id = register_id(child, ids_found, ids_found[p1_id]['level'] + 1, adding_child=True)
                child_ids.append(id)
            ids_found[p1_id]['children'] = child_ids
            ids_found[p2_id]['children'] = child_ids
            return parent1, parent2, children


        for i, line in enumerate(text.split('|')):
            level = i if i == 0 else None
            ret = parse_pedigree_text_line(line, ids_found, level)

        j = convert_to_json_format(ids_found)
        p, act_gens, consistent_per_inferrer, possible_genotypes_per_inferrer = Pedigree.pedigree_from_json(j)
        return j, p, act_gens, consistent_per_inferrer, possible_genotypes_per_inferrer

    def initialise_genome(self):
        g = Gene(AlleleSet.default_alleleset_from_symbol(self.symbol), 1000000,
                 inheritance_pattern=self.inh_type)

        gene_list = [g]
        c = ChromosomeTemplate('PedAutoChrom', 2000000,
                               gene_list if self.chrom_type == ChromosomeTemplate.AUTOSOMAL else [])
        c_X = ChromosomeTemplate('PedXChrom', 2000000, gene_list if self.chrom_type == ChromosomeTemplate.X else [],
                                 type=ChromosomeTemplate.X)
        c_Y = ChromosomeTemplate('PedYChrom', 2000000, gene_list if self.chrom_type == ChromosomeTemplate.Y else [],
                                 type=ChromosomeTemplate.Y)

        return GenomeTemplate(ploidy=2, chromosome_templates=[c], X_chromosome_template=c_X,
                                 Y_chromosome_template=c_Y, name='PedGT')

    def generate(self, hom_rec_partners = False):

        # if self.inh_type != Gene.INH_PATTERN_RECESSIVE:
        #     raise Exception('pedigree type not allowed: ' + self.inh_type)

        #print(str(self.gt))

        adam = self.add_organism(sex=Genome.MALE)
        # eve = Organism.organism_with_random_genotype(gt, sex=Genome.FEMALE, counter_id=next_id)
        # print(str(eve))
        # next_id += 1

        for i in range(0, self.max_levels - 1):
            orgs = self.all_orgs_in_pedigree(include_inlaws=False, level=i)
            for org in orgs:
                r = random.randint(0, 1) # prob of mating
                if i == 0:
                    min_children = 2
                else:
                    min_children = 0
                num_ch = random.randint(min_children, self.max_children)

                if (i == 0 or r <= self.prob_mate)\
                        and (num_ch > 0): # N.B. don't bother mating if no children

                    org_mate = self.add_organism(sex = Genome.MALE if org.sex == Genome.FEMALE else Genome.FEMALE, hom_rec=hom_rec_partners)
                    org.set_partner(org_mate)
                    org.mate(times=num_ch, next_id=self.next_id)
                    self._organisms.extend(org.children)
                    self.next_id += org.num_children


        return adam


class GenotypeInferrer:
    def __init__(self, pedigree, chrom_type, inh_type, def_male_unaf=None, def_female_unaf=None,def_male_af=None, def_female_af=None):
        self.pedigree = pedigree
        self.chrom_type = chrom_type
        self.inh_type = inh_type
        self.def_male_unaf = def_male_unaf
        self.def_female_unaf = def_female_unaf
        self.def_male_af = def_male_af
        self.def_female_af = def_female_af

        self.all_possible_genotypes = self.init_all_possible_genotypes()

    def infer_org(self):

        return False


    @property
    def inferrer_type(self):
        return self.chrom_type + str(self.inh_type)


    def init_all_possible_genotypes(self):
        return []

    def initialise(self):

        orgs = self.pedigree.all_orgs_in_pedigree()
        for org in orgs:
            #affl = org.genome.phenotype_afflicted()[0]
            if org.sex == Genome.MALE:
               if self.is_afflicted(org): # affl.islower():
                   org.set_possible_genotype(self.inh_type, self.chrom_type, self.def_male_af)
               else:
                   org.set_possible_genotype(self.inh_type, self.chrom_type, self.def_male_unaf)
            else:
                if self.is_afflicted(org): #affl.islower():
                    org.set_possible_genotype(self.inh_type, self.chrom_type, self.def_female_af)
                else:
                    org.set_possible_genotype(self.inh_type, self.chrom_type, self.def_female_unaf)



    def is_afflicted(self, org):
        affl_list = org.genome.phenotype_afflicted()
        if len(affl_list) == 0:
            #TODO - replace with actual allele, not hard code A
            affl = 'A'
        else:
            affl = affl_list[0]

        return affl.islower()

    def infer(self):

        self.initialise()

        orgs = self.pedigree.all_orgs_in_pedigree()
        done = False

        consistent = True
        err_msg = None

        while not done:
            num_inferred_this_pass = 0
            for org in orgs:
                changed_this_pass, err_msg_this_pass = self.infer_org(org)
                if changed_this_pass:
                    num_inferred_this_pass += 1
                if err_msg_this_pass is not None:
                    consistent = False
                    err_msg = err_msg_this_pass
                    done = True
                    break
            if num_inferred_this_pass == 0:
                done = True  # Can't do any more

        return consistent, err_msg

    def set_err_msg(self, rule, org):
        return {'rule': rule, 'org': org.counter_id}

    def num_parents_afflicted(self, org):
        _num_parents_afflicted = 0
        for parent in org.parents:
            if self.is_afflicted(parent):# parent.genome.phenotype_afflicted()[0].islower():
                _num_parents_afflicted += 1
        return _num_parents_afflicted

    def show_inferred(self):
        orgs = self.pedigree.all_orgs_in_pedigree()
        inferred_list = []
        key = self.chrom_type + str(self.inh_type)
        for org in orgs:
            if key in org._possible_genotypes:
                inferred = org._possible_genotypes[key]
            else:
                inferred = None
            inferred_list.append(str(org.counter_id) + ':' + ('*' if inferred is None else inferred))
        print(inferred_list)


    def rule_1(self, org):
        return False, False

    def rule_2(self, org):
        return False, False

    def rule_3(self, org):
        return False, False

    def rule_4(self, org):
        return False, False

    def rule_5(self, org):
        return False, False

    def rule_6(self, org):
        return False, False

    def rule_7(self, org):
        return False, False

    def rule_8(self, org):
        return False, False

    def rule_9(self, org):
        return False, False

    def rule_10(self, org):
        return False, False

    def rule_11(self, org):
        return False, False

    def rule_12(self, org):
        return False, False

    def rule_13(self, org):
        return False, False

    def rule_14(self, org):
        return False, False

    def rule_15(self, org):
        return False, False

    def rule_16(self, org):
        return False, False

    def infer_org(self, org):

        any_change = False
        err_msg = None

        # Rule 1  MCA  / MPA
        if org.sex == Genome.MALE and self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_1(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(1, org)

        # Rule 2  MCA  / FPA
        if org.sex == Genome.MALE and self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_2(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(2, org)

        # Rule 3  MCA  / BPA
        if org.sex == Genome.MALE and self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_3(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(3, org)

        # Rule 4  MCA  / NPA
        if org.sex == Genome.MALE and self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_4(org)

            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(4, org)


        # Rule 5  FCA  / MPA
        if org.sex == Genome.FEMALE and self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_5(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(5, org)

        # Rule 6  FCA  / FPA
        if org.sex == Genome.FEMALE and self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_6(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(6, org)

        # Rule 7  FCA  / BPA
        if org.sex == Genome.FEMALE and self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_7(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(7, org)

        # Rule 8  FCA  / NPA
        if org.sex == Genome.FEMALE and self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_8(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(8, org)


        # Rule 9  MCN  / MPA
        if org.sex == Genome.MALE and not self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_9(org)

            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(9, org)

        # Rule 10  MCN  / FPA
        if org.sex == Genome.MALE and not self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_10(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(10, org)

        # Rule 11  MCN  / BPA
        if org.sex == Genome.MALE and not self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):

                changed, err = self.rule_11(org)
                if changed:
                    any_change = True
                if err:
                    err_msg = self.set_err_msg(11, org)

        # Rule 12  MCN  / NPA
        if org.sex == Genome.MALE and not self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):

                changed, err = self.rule_12(org)
                if changed:
                    any_change = True
                if err:
                    err_msg = self.set_err_msg(12, org)

         # Rule 13  FCN  / MPA
        if org.sex == Genome.FEMALE and not self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and not self.is_afflicted(org.female_parent):
            changed, err = self.rule_13(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(13, org)

        # Rule 14  FCN  / FPA
        if org.sex == Genome.FEMALE and not self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_14(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(14, org)

        # Rule 15  FCN  / BPA
        if org.sex == Genome.FEMALE and not self.is_afflicted(org) \
                and org.has_parents \
                and self.is_afflicted(org.male_parent) and self.is_afflicted(org.female_parent):
            changed, err = self.rule_15(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(15, org)

        # Rule 16  FCN  / NPA
        if org.sex == Genome.FEMALE and not self.is_afflicted(org) \
                and org.has_parents \
                and not self.is_afflicted(org.male_parent) and  not self.is_afflicted(org.female_parent):
            changed, err = self.rule_16(org)
            if changed:
                any_change = True
            if err:
                err_msg = self.set_err_msg(16, org)


        return any_change, err_msg



class ARGenotypeInferrer(GenotypeInferrer):

    def __init__(self, pedigree):

        #affl = pedigree.all_orgs_in_pedigree()[0].genome.phenotype_afflicted()[0]
        alleles = pedigree.all_orgs_in_pedigree()[0].afflicted
        if len(alleles) == '':
            # TODO replace with proper sourcing of allele name
            self.allele = 'A'
        else:
            self.allele = alleles[0]
        #allele = pedigree.all_orgs_in_pedigree()[0].genotype()[0]
        def_male_unaf = None
        def_female_unaf = None
        def_male_af = self.allele.lower() + self.allele.lower()
        def_female_af = self.allele.lower() + self.allele.lower()

        super(ARGenotypeInferrer, self).__init__(pedigree,
                                                 ChromosomeTemplate.AUTOSOMAL,
                                                 Gene.INH_PATTERN_RECESSIVE,
                                                 def_male_unaf=def_male_unaf,
                                                 def_female_unaf=def_female_unaf,
                                                 def_male_af=def_male_af,
                                                 def_female_af=def_female_af
                                                 )

    def init_all_possible_genotypes(self):
        return [
            self.allele.upper() + self.allele.lower(),
            self.allele.lower() + self.allele.lower(),
            self.allele.upper() + '-'
        ]


    def rule_1(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_2(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_4(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        if err:
            return changed, err

        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())

        return changed, err

    def rule_5(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_6(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_8(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())

        if err:
            return changed, err

        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                           self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_9(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_10(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_11(self, org):
        return False , True

    def rule_13(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_14(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_15(self, org):
        return False, True



class ADGenotypeInferrer(GenotypeInferrer):

    def __init__(self, pedigree):

        alleles = pedigree.all_orgs_in_pedigree()[0].afflicted
        if len(alleles) == '':
            #TODO replace with proper sourcing of allele name
            self.allele = 'A'
        else:
            self.allele = alleles[0]
            #pedigree.all_orgs_in_pedigree()[0].genotype()[0] #pedigree.all_orgs_in_pedigree()[0].genome.phenotype_afflicted()[0]
        def_male_unaf = self.allele.upper() + self.allele.upper()
        def_female_unaf = self.allele.upper() + self.allele.upper()
        def_male_af = None
        def_female_af = None

        super(ADGenotypeInferrer, self).__init__(pedigree,
                                                 ChromosomeTemplate.AUTOSOMAL,
                                                 Gene.INH_PATTERN_DOMINANT,
                                                 def_male_unaf=def_male_unaf,
                                                 def_female_unaf=def_female_unaf,
                                                 def_male_af=def_male_af,
                                                 def_female_af=def_female_af
                                                 )

    def init_all_possible_genotypes(self):
        return [
            self.allele.upper() + self.allele.upper(),
            self.allele.upper() + self.allele.lower(),
            '-' + self.allele.lower()
        ]

    def rule_1(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_2(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_4(self, org):
        return False, True

    def rule_5(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_6(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_8(self, org):
        return False, True

    def rule_9(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_10(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_11(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        if err:
            return changed, err

        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())

        return changed, err

    def rule_13(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_14(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())
        return changed, err

    def rule_15(self, org):
        changed, err = org.male_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                             self.allele.upper() + self.allele.lower())
        if err:
            return changed, err


        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               self.allele.upper() + self.allele.lower())
        return changed, err



class XRGenotypeInferrer(GenotypeInferrer):

    def __init__(self, pedigree):

        #affl = pedigree.all_orgs_in_pedigree()[0].genome.phenotype_afflicted()[0]
        alleles = pedigree.all_orgs_in_pedigree()[0].afflicted
        if len(alleles) == '':
            # TODO replace with proper sourcing of allele name
            self.allele = 'A'
        else:
            self.allele = alleles[0]
        #allele = pedigree.all_orgs_in_pedigree()[0].genotype()[0]
        def_male_unaf = 'X' + self.allele.upper() + 'Y'
        def_female_unaf = None
        def_male_af = 'X' + self.allele.lower() + 'Y'
        def_female_af = 'X' + self.allele.lower() + 'X' + self.allele.lower()

        super(XRGenotypeInferrer, self).__init__(pedigree,
                                                 ChromosomeTemplate.X,
                                                 Gene.INH_PATTERN_RECESSIVE,
                                                 def_male_unaf=def_male_unaf,
                                                 def_female_unaf=def_female_unaf,
                                                 def_male_af=def_male_af,
                                                 def_female_af=def_female_af
                                                 )

    def init_all_possible_genotypes(self):
        return [
            'X' + self.allele.upper() + 'X' + self.allele.lower(),
            'X' + self.allele.upper() + 'X' + '-',
            'X' + self.allele.lower() + 'X' + self.allele.lower(),
            'X' + self.allele.upper() + 'Y',
            'X' + self.allele.lower() + 'Y'
        ]

    def rule_1(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err


    def rule_4(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               'X' + self.allele.upper() + 'X' + self.allele.lower())

        return changed, err

    def rule_5(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err

    def rule_6(self, org):
       return False, True

    def rule_8(self, org):
        return False, True


    def rule_10(self, org):
        return False, True

    def rule_11(self, org):
        return False, True

    def rule_13(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err

    def rule_14(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err

    def rule_15(self, org):
        return False, True


class XDGenotypeInferrer(GenotypeInferrer):

    def __init__(self, pedigree):

        #affl = pedigree.all_orgs_in_pedigree()[0].genome.phenotype_afflicted()[0]
        alleles = pedigree.all_orgs_in_pedigree()[0].afflicted
        if len(alleles) == '':
            # TODO replace with proper sourcing of allele name
            self.allele = 'A'
        else:
            self.allele = alleles[0]
        #allele = pedigree.all_orgs_in_pedigree()[0].genotype()[0]
        def_male_unaf = 'X' + self.allele.upper() + 'Y'
        def_female_unaf = 'X' + self.allele.upper() + 'X' + self.allele.upper()
        def_male_af = 'X' + self.allele.lower() + 'Y'
        def_female_af = None

        super(XDGenotypeInferrer, self).__init__(pedigree,
                                                 ChromosomeTemplate.X,
                                                 Gene.INH_PATTERN_DOMINANT,
                                                 def_male_unaf=def_male_unaf,
                                                 def_female_unaf=def_female_unaf,
                                                 def_male_af=def_male_af,
                                                 def_female_af=def_female_af
                                                 )

    def init_all_possible_genotypes(self):
        return [
            'X' + self.allele.upper() + 'X' + self.allele.upper(),
            'X' + '-' + 'X' + self.allele.lower(),
            'X' + self.allele.upper() + 'X' + self.allele.lower(),
            'X' + self.allele.upper() + 'Y',
            'X' + self.allele.lower() + 'Y'
        ]

    def rule_1(self, org):
        return False, True

    def rule_4(self, org):
        return False, True

    def rule_5(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err

    def rule_6(self, org):
        changed, err = org.set_possible_genotype(self.inh_type, self.chrom_type,
                                                 'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err

    def rule_8(self, org):
        return False, True

    def rule_10(self, org):
        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               'X' + self.allele.upper() + 'X' + self.allele.lower())
        return changed, err


    def rule_11(self, org):

        changed, err = org.female_parent.set_possible_genotype(self.inh_type, self.chrom_type,
                                                               'X' + self.allele.upper() + 'X' + self.allele.lower())

        return changed, err

    def rule_13(self, org):
        return False, True

    def rule_15(self, org):
        return False, True


class YGenotypeInferrer(GenotypeInferrer):

    def __init__(self, pedigree):

        #affl = pedigree.all_orgs_in_pedigree()[0].genome.phenotype_afflicted()[0]
        alleles = pedigree.all_orgs_in_pedigree()[0].afflicted
        if len(alleles) == '':
            # TODO replace with proper sourcing of allele name
            self.allele = 'A'
        else:
            self.allele = alleles[0]
        #allele = pedigree.all_orgs_in_pedigree()[0].genotype()[0]
        def_male_unaf =  'XY' # 'Y' + self.allele.upper() + 'X'
        def_female_unaf = 'XX' #'X'  + 'X'
        def_male_af =  'Xy'  #'Y' + self.allele.lower() + 'X'
        def_female_af = None #Not allowed

        super(YGenotypeInferrer, self).__init__(pedigree,
                                                 ChromosomeTemplate.Y,
                                                 Gene.INH_PATTERN_RECESSIVE,
                                                 def_male_unaf=def_male_unaf,
                                                 def_female_unaf=def_female_unaf,
                                                 def_male_af=def_male_af,
                                                 def_female_af=def_female_af
                                                 )

    def init_all_possible_genotypes(self):
        return [
            'X' + 'X',
            'X' + 'Y',
            'X' + 'y'
        ]

    def infer_org(self, org):

        any_change = False
        err_msg = None

        # Rule 0 FCA
        if org.sex == Genome.FEMALE and self.is_afflicted(org):
            err_msg = self.set_err_msg(0, org)

        # Rule 17 Male afflicted / male parent unafflicted or Male unafflicted / male parent afflicted
        if   (org.sex == Genome.MALE and self.is_afflicted(org) \
              and  org.has_parents \
              and not self.is_afflicted(org.male_parent)) or \
             (org.sex == Genome.MALE and not self.is_afflicted(org) \
              and  org.has_parents \
              and self.is_afflicted(org.male_parent)):
             err_msg = self.set_err_msg(17, org)


        return any_change, err_msg


if __name__ == '__main__':
    debug = 0

    random.seed(42)
    c1 = ChromosomeTemplate.from_symbol_list(['A-R-1000', 'B-D-2000', 'C-R-100000000', 'D-R-200000000'], 'Chrom 1')
    gt = GenomeTemplate(ploidy=2,chromosome_templates = [c1], name='Onechrom')
    tst_1 = Organism.organism_with_het_genotype(gt)
    tst_2 = Organism.organism_with_het_genotype(gt)
    children = tst_1.mate(tst_2, times=4)


    ret = Pedigree.pedigree_from_text('1MU-2FU:3MU,4FU,5FU|3-6FA:7FU,8MA,9MA|9-13FA:14MA,15MA,16FA|5FU-10MU:11FU,12FU|11-17MA:18FU,19FU,20FU,21FU')

    ret = Pedigree.pedigree_from_text('1M-2F:3F,4MA|3-5MA:6F,7MA,8F')

    from scipy.stats import chisquare
    # obs = [1469, 138, 5]
    # exp = [1467.4, 141.2, 3.4]
    obs = [1489, 108, 15]
    exp = [1467.4, 141.2, 3.4]
    chisq, p = chisquare(obs, exp, ddof=1)
    F = 1 - (obs[1] / exp[1])

    phens = Genome.test_cross_het_gametes_to_phenotypes(['ABC', 'ABc'])

    org = Organism.organism_from_gametes({'ABC':37,'ABc':378,'AbC':10, 'Abc':100,'aBC':88,'aBc':7, 'abC': 344,'abc':36})
    #a1 = Allele('A')
    #a2 = Allele('a')
    alleles_a = AlleleSet.default_alleleset_from_symbol('A')
    g1 = Gene(alleles_a,10000000)
    g1_dom = Gene(alleles_a,10000000,inheritance_pattern=Gene.INH_PATTERN_DOMINANT)

    g2 = Gene(AlleleSet.default_alleleset_from_symbol('B'),170)
    g3 = Gene(AlleleSet.default_alleleset_from_symbol('C'),80000000)
    g5 = Gene(AlleleSet.default_alleleset_from_symbol('D'),50000000)
    g4 = Gene(AlleleSet.default_alleleset_from_symbol('G'),150)
    g6 = Gene(AlleleSet.default_alleleset_from_symbol('E'),30000000)
    g6_dom = Gene(AlleleSet.default_alleleset_from_symbol('E'),30000000, inheritance_pattern=Gene.INH_PATTERN_DOMINANT)

    g7 = Gene(AlleleSet.default_alleleset_from_symbol('F'),90000000)
    g7_dom = Gene(AlleleSet.default_alleleset_from_symbol('F'), 90000000,inheritance_pattern = Gene.INH_PATTERN_DOMINANT )

    print(str(alleles_a))
    print(g1)

    d = g1._to_attr_dict()
    g_new = Gene._from_attr_dict(d)

    d = g3._to_attr_dict()

    g_new = Gene._from_attr_dict(d)

    c1 = ChromosomeTemplate('3',200,[g2])

    c2 = ChromosomeTemplate('XL',350,[g1,g3,g5])

    c_ped = ChromosomeTemplate('XL', 1000, [g1])

    c_attr_dict = c2._to_attr_dict()

    c2_inflated = ChromosomeTemplate._from_attr_dict(c_attr_dict)


    c3 = ChromosomeTemplate('XR',1000,[g4])
    print(c1)

    c4 = ChromosomeTemplate('XChrom',30000000,[g6],type=ChromosomeTemplate.X)
    c5 = ChromosomeTemplate('YChrom', 10000000, [g7], type=ChromosomeTemplate.Y)


    gt = GenomeTemplate(ploidy=2,chromosome_templates = [c1,c2,c3], X_chromosome_template=c4, Y_chromosome_template=c5, name='Fivechroms')
    print(str(gt))

    c_ped_x = ChromosomeTemplate('XChrom',30000000,[g6],type=ChromosomeTemplate.X)
    c_ped_y = ChromosomeTemplate('YChrom', 10000000, [g7], type=ChromosomeTemplate.Y)
    gt_ped = GenomeTemplate(chromosome_templates = [c_ped],X_chromosome_template = c_ped_x, Y_chromosome_template = c_ped_y, name='Pedigree' )
    # org_ped = Organism.organism_with_random_genotype(gt_ped)
    # print('org_ped: ',org_ped)
    # org_ped_2 = Organism.organism_with_random_genotype(gt_ped)
    # print('org_ped_2: ',org_ped_2
    # org_ped.set_partner(org_ped_2)
    #
    # org_ped.mate(times=4)
    # print('org_ped child 1: ',org_ped.children[0])
    c_ped_dom = ChromosomeTemplate('XL', 1000, [g1_dom])
    c_ped_x_dom = ChromosomeTemplate('XChrom', 30000000, [g6_dom], type=ChromosomeTemplate.X)
    c_ped_y_dom = ChromosomeTemplate('YChrom', 10000000, [g7_dom], type=ChromosomeTemplate.Y)

    gt_tst = GenomeTemplate(chromosome_templates=[c_ped_dom], X_chromosome_template=c_ped_x_dom, Y_chromosome_template=c_ped_y_dom,
                            name='Pedigree')

    tst = Organism.organism_with_het_genotype(gt_tst,sex=Genome.MALE,counter_id=-3)

    next_id = 1
    org_ped_root_1 = Organism.organism_with_hom_dominant_genotype(gt_ped, sex=Genome.MALE, counter_id=next_id)
    print(str(org_ped_root_1))
    next_id +=1
    print(org_ped_root_1.genome.phenotype_afflicted())

    org_ped_root_2 =Organism.organism_with_het_genotype(gt_ped, sex=Genome.FEMALE, counter_id=next_id)
    org_ped_root_1.set_partner(org_ped_root_2)
    next_id +=1
    org_ped_root_1.mate(times=4, next_id=next_id)
    next_id += org_ped_root_2.num_children

    org_ped_1_1 = org_ped_root_1.children[0]
    print(org_ped_1_1.genome.phenotype_afflicted())
    print(org_ped_1_1.genome.sex_pair)
    org_ped_mate = Organism.organism_with_hom_recessive_genotype(gt_ped, sex=Genome.MALE if org_ped_1_1.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_1_1.set_partner(org_ped_mate)
    org_ped_1_1.mate(times=2, next_id=next_id)
    next_id +=org_ped_1_1.num_children

    org_ped_2_1  = org_ped_1_1.children[0]
    org_ped_mate = Organism.organism_with_hom_recessive_genotype(gt_ped, sex=Genome.MALE if org_ped_2_1.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_2_1.set_partner(org_ped_mate)

    org_ped_1_2 = org_ped_root_1.children[1]
    org_ped_mate = Organism.organism_with_het_genotype(gt_ped, sex=Genome.MALE if org_ped_1_2.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_1_2.set_partner(org_ped_mate)
    org_ped_1_2.mate(times=1, next_id=next_id)
    next_id +=org_ped_1_2.num_children

    org_ped_2_2  = org_ped_1_2.children[0]
    org_ped_mate = Organism.organism_with_hom_recessive_genotype(gt_ped, sex=Genome.MALE if org_ped_2_2.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_2_2.set_partner(org_ped_mate)
    org_ped_2_2.mate(times=2, next_id=next_id)
    next_id +=org_ped_2_2.num_children


    org_ped_1_3 = org_ped_root_1.children[2]
    org_ped_1_3_mate = Organism.organism_with_hom_recessive_genotype(gt_ped, sex=Genome.MALE if org_ped_1_3.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_1_3.set_partner(org_ped_1_3_mate)
    org_ped_1_3.mate(times=1, next_id = next_id)
    next_id += org_ped_1_3.num_children

    org_ped_2_3 = org_ped_1_3.children[0]
    org_ped_2_3_mate = Organism.organism_with_hom_recessive_genotype(gt_ped, sex=Genome.MALE if org_ped_2_3.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_2_3.set_partner(org_ped_2_3_mate)
    org_ped_2_3.mate(times=4, next_id = next_id)
    next_id += org_ped_2_3.num_children


    org_ped_1_4 = org_ped_root_1.children[3]
    org_ped_1_4_mate = Organism.organism_with_hom_dominant_genotype(gt_ped, sex=Genome.MALE if org_ped_1_4.sex == Genome.FEMALE else Genome.FEMALE, counter_id=next_id)
    next_id +=1
    org_ped_1_4.set_partner(org_ped_1_4_mate)
    org_ped_1_4.mate(times=4, next_id = next_id)
    next_id += org_ped_1_4.num_children

    print(org_ped_1_4.descendants)
    gt_attr_dict = gt._to_attr_dict()

    print('Descendant tree')
    print('\t' * org_ped_root_1.level + str(org_ped_root_1) + ' affl: ' + str(org_ped_root_1.genome.phenotype_afflicted()))
    org_ped_root_1.print_org_tree_below()

    print(org_ped_root_1.orgs_afflicted_below('a'))
    print(org_ped_root_1.orgs_afflicted_below('e'))
    print(org_ped_root_1.orgs_afflicted_below('f'))

#    x= org_ped_1_4.all_orgs_in_pedigree(include_in_laws =True, level=4)

    #adam = Organism.generate_pedigree(max_levels=4)

    prob_mate = 0.9 #0.5
    max_children = 8 #4
    max_levels = 5 # 4
    p = Pedigree(max_levels=max_levels, inh_type=Gene.INH_PATTERN_RECESSIVE, chrom_type=ChromosomeTemplate.X, prob_mate = prob_mate, max_children = max_children)
    adam = p.generate(hom_rec_partners=False)

    j = p.to_json()
    next_l = adam.next_older_sibling

    neph = adam._get_nephews()

    neph_2 = adam.children[0]._get_nephews()

    adam.print_org_tree_below()

    inferrers = [ARGenotypeInferrer(p), ADGenotypeInferrer(p), XRGenotypeInferrer(p), XDGenotypeInferrer(p), YGenotypeInferrer(p)]
    for inferrer in inferrers:
        consistent, err_msg = inferrer.infer()
        if consistent:
            print(f'Consistent with {inferrer.inferrer_type}')
        else:
            print(f'Inconsistent with {inferrer.inferrer_type} -  {err_msg}')
        inferrer.show_inferred()

    act_gens = []
    for org in p.all_orgs_in_pedigree():
        act_gens.append(f'{org}')
    print('Act: ')
    print(str(act_gens))


    consistent_per_inferrer = {}
    for i in range(50):
        print('i: ',i)
        p = Pedigree(max_levels=max_levels, inh_type=Gene.INH_PATTERN_RECESSIVE, chrom_type=ChromosomeTemplate.Y, prob_mate = prob_mate, max_children = max_children)
        adam = p.generate(hom_rec_partners=False)

        inferrers = [ARGenotypeInferrer(p), ADGenotypeInferrer(p), XRGenotypeInferrer(p), XDGenotypeInferrer(p), YGenotypeInferrer(p)]
        for inferrer in inferrers:
            consistent, err_msg = inferrer.infer()

            if consistent:
                if inferrer.inferrer_type in consistent_per_inferrer:
                    consistent_per_inferrer[inferrer.inferrer_type] +=1
                else:
                    consistent_per_inferrer[inferrer.inferrer_type] = 1
            #     print(f'Consistent with {inferrer.inferrer_type}')



    gt_prob_test = GenomeTemplate(ploidy=2,chromosome_templates = [c2], name='Probtest')
    org_prob_test = Organism.organism_with_het_genotype(gt_prob_test, rand_phase=True)
    print(org_prob_test.genome.possible_gametes_formatted())
    print('org prob test genotype: ',org_prob_test.genotype())


    org = Organism.organism_with_random_genotype(gt)
    print(org)

    org_attr_dict = org._to_attr_dict()


    ch1 = ChromosomeTemplate('3',100,[g1])
    ch2 = ChromosomeTemplate('XL-group1',300,[g2])
    gt2 = GenomeTemplate(ploidy=2,chromosome_templates = [ch1,ch2],name='Twochroms')
    gt2_attr_dict = gt2._to_attr_dict()
    gt2_inflated = GenomeTemplate._from_attr_dict(gt2_attr_dict)
    print(gt2_inflated)

    org1 = Organism.organism_with_random_genotype(gt2)
    print('org1: ',org1)

    org1_deflated = org1._to_attr_dict()
    org1_inflated = Organism._from_attr_dict(org1_deflated)
    print('org1 inflated: ',org1_inflated)

    print('gam1: ',org1.genome.get_parental_gamete(0))
    print('gam2: ', org1.genome.get_parental_gamete(1))

    org2 = Organism.organism_with_random_genotype(gt2)
    print('org2: ',org2)
    child = org1.mate(org2)
    print('child: ',child)


    print('org1: \n',org1.genotype())
    print('org2: \n',org2.genotype())

    new = org1.mate(org2)
    print('new:\n', new.genotype())

    org11 = Organism.organism_with_random_genotype(gt, sex=Genome.MALE)
    org22 = Organism.organism_with_random_genotype(gt, sex=Genome.FEMALE)

    org11.genome.chromosome_pairs[1].possible_gametes()
    print(str(org22))
    print('Org11')
    print(org11.genotype())
    print('Org22')
    print(org22.genotype())
    new33 = org11.mate(org22)
    print('New33')
    print(new33.genotype())
    print(new33)

    children = []
    for i in range(50):
        children.append(org11.mate(org22))

    print('A\u00B2') #Superscript
    print('A\u2082') #subscript
    print('X\u00BFX\u208f')
    print('org11: ',org11)
    print('org22: ',org22)
    if debug > 0:
        print('Children:')
        for child in children:
         print(child)

    gt_linkage = GenomeTemplate(ploidy=2, chromosome_templates=[c2])
    print(gt)

    org_hom = Organism.organism_with_hom_recessive_genotype(gt_linkage)
    print(org_hom)
    print (org_hom.genotype())

    org_het = Organism.organism_with_het_genotype(gt_linkage, rand_phase=True)
    print(org_het)
    print (org_het.genotype())

    print('het gam1: ', org_het.genome.get_parental_gamete(0))
    print('het gam2: ', org_het.genome.get_parental_gamete(1))

    children = []
    for i in range(1000):
        children.append(org_het.mate(org_hom))

    print('child1 gam1: ', children[0].genome.get_parental_gamete(0,sort_alpha=True))
    print('child1 gam2: ', children[0].genome.get_parental_gamete(1,sort_alpha=True))
    print('child2 gam1: ', children[1].genome.get_parental_gamete(0,sort_alpha=True))
    print('child2 gam2: ', children[1].genome.get_parental_gamete(1,sort_alpha=True))
    genotypes = {}
    for child in children:
        genotype = child.genotype()
        if genotype in genotypes:
           genotypes[genotype] +=1
        else:
            genotypes[genotype] = 1
        #print(child)

    print('org het: ',org_het.genotype())
    print('org hom: ',org_hom.genotype())
    print(genotypes)

    print (org_hom.genome.chromosome_pairs[0].phenotype())