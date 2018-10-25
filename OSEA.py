import pandas as pd
import numpy as np
import stats_test 
import biom

class OSEA(object):
    """ OTU Set Enrichment Analysis,like GSEA.
    Here,Set are chose by user,it can be phylum genus etc.

    Attributes:

    """

    def __init__(self,rank_list, Taxon_file='Taxon.tsv', set_level='Phylum'):    
        """ Init OSEA class
        Args:
            rank_list: It is a dictionary which store the some values(corr) in 
            order. e.g.
                {'otu1':0.9,'otu2':0.8,'otu3:0.2}
            taxon: include the lineage information about OTUs.
            set_level: it can be species ,genus, family, phylum etc.
        """
        self.rank_list = rank_list
        self.taxon = pd.read_csv(Taxon_file, sep='\t')
        self.set_level = set_level
        self.sets = self.get_sets()
        self.get_ES()

    def get_sets(self):
        """ get set according to self.set_level
        Args:
        Return: a list of set.For example:

        ['phylum0','phylum1','phylum2']
        """
        mapped_level={'level0':'Kingdom','level1':'Phylum','level2':'Class',
                      'level3':'Order','level4':'Family','level5':'Genus',
                      'level6':'Species','level7':'OTU'}
        index = 1
        for level,name in mapped_level.items():
            if  name == self.set_level:
                index = int(level[-1])
                self.index = index
        sets = []  # add the phylum name to this sets.
        for ele in self.taxon['Taxon']:
            tmp = ele.split(';')
            if len(tmp) > index+1:
                unit = tmp[index]
                if unit in sets:
                    pass
                else:
                    sets.append(unit)
        return sets  # return 
                
        
    def get_ES(self, power=1):
        """ Compute the enrichment score for every set.
        Arg:
            power: power number used for compute es in the iterator step.
        """
        es = {}  # enrichment score
        rank_list = self.rank_list
        index = self.index
        for ele in self.sets:  # OTU set
            tmp = 0
            tmp_unit = ''
            p_hit = [0] # enrichment score for every score.p_hit minus p_miss
            p_miss = [0] 
            rank_pow_sum = 0
            for otu in rank_list:
                rank_pow_sum += rank_list[otu]**power
            i = 1 # index used for rank_list
            N = len(rank_list)
            for otu in rank_list:  # go through rank list
                taxo = self.taxon.loc[self.taxon['Feature ID']==otu]['Taxon']
                for elem in taxo:
                    f = elem.split(';')
                    if len(f) > index:
                        tmp_unit = f[index]
                    else:
                        tmp_unit = ''
                N_h = 0
                if tmp_unit == ele:  # hit, like gene in some Set
                    N_h += 1
                    assert rank_pow_sum != 0
                    tmp_p_hit = p_hit[-1]+(rank_list[otu]**power)/rank_pow_sum
                    p_hit.append(tmp_p_hit)
                    tmp_p_miss = (i-N_h)/(N-N_h)  # TODO
                    p_miss.append(tmp_p_miss)
                else:  # miss
                    tmp_p_hit = p_hit[-1]
                    p_hit.append(tmp_p_hit)
                    tmp_p_miss = (i-N_h)/(N-N_h)  # TODO
                    p_miss.append(tmp_p_miss)
                i += 1
            assert len(p_hit) == len(p_miss)
            point_es = []
            for j in range(len(p_hit)):
                point_es.append(p_hit[j]-p_miss[j])
            self.point_es = point_es
            m = max(point_es)
            absm = max(map(abs,point_es))
            if absm > m:
                es = -absm
            else:
                es = m
        self.es = es

def permutation_to_obtain_ranklist(feature_table, test_method_name='t_test'):
    """
    Randomly generate 0,1 labels , rank list 
    Arg:
        feature_table: it is the file name of the biom format otu table.
        test_method_name: t_test or F_test  a string.
    Return: 
        rank_list:a p value dict.
    """
    df = biom.load_table(feature_table).to_dataframe().transpose().to_dense()
    n = df.shape[0]
    part1_index = []
    part2_index = []
    for i in range(n):
        if np.random.randint(0,2):
            part1_index.append(i)
        else:
            part2_index.append(i)
    part1 = df.iloc[part1_index, :] 
    part2 = df.iloc[part2_index, :] 
    methods = {'t_test': stats_test.t_test,
               'F_test': stats_test.F_test}
    test_method = methods[test_method_name]
    rank_list_unsort = {}
    for col in part1.columns:
        tmp_pvalue = test_method(part1[col],part2[col])
        rank_list_unsort[col]=tmp_pvalue
    tmp_list = sorted(rank_list_unsort, key=rank_list_unsort.get,reverse=True)
    rank_list = {}
    for ele in tmp_list:
        rank_list[ele]=rank_list_unsort[ele]
    return rank_list
def find_index(index,indexes):
    """find the index between indexes.
    For example,indexes = [1,2,3],index = 1.3
    it will return the index 0
    """
    for i in range(len(indexes)):
        if i == len(indexes)-1:
            return i
        left = index -indexes[i]
        right = indexes[i+1]-index
        if left <= 0:
            return 0
        if left > 0 and right > 0:
            if left <= right:
                return i
            else:
                return i+1
def p_value(sample_value,distribution):
    """Obtain the pvalue when given the sample value and the distribution
    Args:
        sample_value:sample_value(a numerical value)
        distribution: a dict saving values and their probability.
    Return:
        pvalue
    """
    index = find_index(sample_value,distribution[1])
    prob = distribution[0][index]
    pvalue = 0
    for ele in distribution[0]:
        if ele <= prob:
            pvalue += ele
    return pvalue
def generate_distribution(arr):
    """Generate a distribution of the given array.
    Return : two arrays.
        arrays[0] is the probability,arrays[1] is the indexes
        For example,
            [0.1,0.1,0.8],[1,2,3]
    
    """
    values, indexes = np.histogram(arr,bins=3)
    values = values / sum(values)
    return values, indexes
        







