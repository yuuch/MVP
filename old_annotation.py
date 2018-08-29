import biom
import pandas as pd 
from Bio import Phylo 
import plotly
import plotly.graph_objs as go
'''
annotation for a sample
'''
class Annotation(object):
    def __init__(self, tree, feature_table_file, taxo_file):
        self.tree = tree
        print()
        self.feature_table = biom.load_table(feature_table_file).to_dataframe().to_dense()
        self.taxonomy = pd.read_csv(taxo_file,sep='\t')
        self.taxonomy.index = self.taxonomy['Feature ID']
        self.map_taxo(self.tree)

    def map_taxo(self,tree):
        """
        map feature id to taxonomy and compute the number of every phylum,class and so on .
        """
        terminals = tree.get_terminals()
        feature_names = [node.name for node in terminals]
        feature_name_count_dict = {}
        for feature_name in feature_names:
            # feature_count: the number of seq in one feature(otu)
            feature_count = sum(list(self.feature_table.loc[feature_name]))
            if feature_count == 0:
                continue
            taxo_of_feature = self.taxonomy.loc[feature_name]['Taxon'].split(';')
            levels = {}
            for i in range(8):
                levels['level'+str(i)]={}
            #levels['level7'][feature_name] = [feature_count]
            for i in range(len(taxo_of_feature)):
                if taxo_of_feature[i] in levels['level'+str(i)]:
                    levels['level'+str(i)][taxo_of_feature[i]] += feature_count
                else:
                    levels['level'+str(i)][taxo_of_feature[i]] = feature_count
            feature_name_count_dict[feature_name]=levels
        sum_dict = {}
        for i in range(8):
            sum_dict['level'+str(i)]={}
        for feature in feature_name_count_dict:
            for level in feature_name_count_dict[feature]:
                for ele in feature_name_count_dict[feature][level]:
                    if ele in sum_dict[level]:
                        sum_dict[level][ele] += feature_name_count_dict[feature][level][ele]
                    else:
                        sum_dict[level][ele] = feature_name_count_dict[feature][level][ele]
        self.barplot_dict = sum_dict

    """

    performed other place.

    def add_node_num(self,tree):
        nodes = tree.get_terminals()+tree.get_non_terminals()
        i = 0
        for node in nodes:
            node.node_num = i 
            i += 1
    """

    def plot_annotation(self):
        mapped_level={'level0':'Kingdom','level1':'Phylum','level2':'Class','level3':'Order',
                 'level4':'Family','level5':'Genus','level6':'Species','level7':'OTU'}
        data = []
        for level in self.barplot_dict:
            for ele in self.barplot_dict[level]:
                trace = go.Bar(
                        y = [mapped_level[level]],
                        x = [self.barplot_dict[level][ele]],
                        name = ele,
                        orientation = 'h'
                    )
                data.append(trace)
        layout = go.Layout(
            #autosize =False,
            #height = 500,
            #width = 1000,
            barmode = 'stack'
        )
        fig = go.Figure(data=data,layout=layout)
        div_str = plotly.offline.plot(fig,filename='stack_bar.html',output_type='div')
        return div_str
                




        



        
        


