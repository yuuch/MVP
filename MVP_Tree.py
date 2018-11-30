#!/usr/bin/python3
import pandas as pd 
import numpy as np 
import biom 
from Bio import Phylo
class MvpTree(object):
    """ A tree which contain the information from feature table, taxonomy file
    and metadata.
    """
    def __init__(self, feature_table, tree, taxonomy):
        """ Init the class with some basic files generated from qiime2."""
        self.tree = Phylo.read(tree, 'newick')
        self.feature_table = biom.load_table(feature_table).to_dataframe()
        tmp_taxo = pd.read_csv(taxonomy, sep='\t')
        self.taxonomy = tmp_taxo.set_index(tmp_taxo['Feature ID'])
    def generate_sample_dict(self):
        """ For every node in the tree ,obtain a sample dict for it.
        e.g sample_dict = {'S0':11,'S1':12, 'S3': 13}
        here we need the tree and the feature table.
        """
        pass
    def generate_domain_otu(self):
        """We want to obtain the domain node (which posseses the most sequense)
        of every node (used for color in the tree plot)
        """
        pass
    




