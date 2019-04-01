from Bio import Phylo
import biom
import copy 
import math
import numpy as np 
import pandas as pd
import plotly
import re
import scipy.stats as stats

import rectangle_subtree
import stats_test
class MvpTree(object):
    """ a class merged the metadata, otutable and phylogenetic tree.
    """
    def __init__(self, feature_table_path, tree_path, metadata_path, \
        taxonomy_path=None, ID_num=0):
        self.feature_table = biom.load_table(feature_table_path)\
            .to_dataframe()
        self.tree = Phylo.read(tree_path, "newick")
        self.metadata = pd.read_csv(metadata_path, sep='\t')
        self.init_col = self.feature_table.index
        self.get_normalize_feature_table()
        self.metadata_set_index()
        self.get_featured_tree()
        if taxonomy_path:
            self.get_taxonomy_info(taxonomy_path)
            self.get_domain_otu()
        self.get_subtree(ID_num)
        self.get_GI()


    def get_taxonomy_info(self,taxonomy_path):
        """ get taxonomy infomation for every otu."""
        Taxon = 'Taxon'
        p1 = '.*[Tt][Aa][Xx][Oo].*'
        pattern = re.compile(p1)
        feature_id = 'Feature ID'
        p2 = '.*[Ii][Dd].*'
        pattern2 = re.compile(p2)
        
        try:
            taxonomy_df = pd.read_csv(taxonomy_path, sep= '\t')
            print('valid taxonomy file')
        except:
            print('unvalid  taxonomy path')
        for ele in taxonomy_df.columns:
            if len(pattern.findall(ele)) >0:
                if pattern.findall(ele)[0]>3 :
                    Taxon = ele
            else:
                pass
            break
        for ele in taxonomy_df.columns:
            if len(pattern2.findall(ele)) > 0:
                if len(pattern2.findall(ele)[0])>3 :
                    feature_id = ele
            else:
                pass
            break
        taxonomy_df = taxonomy_df.set_index(feature_id)
        self.lineage = taxonomy_df[Taxon]
        
    def get_colors(self, colors, color_index):
        for clade in self.feature_tree.find_clades(order='level'):
            try:
                lineages = self.lineage[clade.domain_otu]
                #print(lineages)
                lineages = lineages.split(';')
                phylumn_name = lineages[1] #phylumn name
                #print (phylumn_name)
                temp_index = color_index[phylumn_name]
                #print(temp_index)
                clade.plot_color = colors[temp_index]
                #print(clade.plot_color)
            except:
                clade.plot_color = 'rgb(0,0,0)'

    def metadata_set_index(self):
        #metadata = self.metadata
        #feature_table = self.feature_table.transpose()
        p1 = '.*[Ss][Aa][Mm][Pp][Ll][Ee].*[IiNn][DdAa]'
        pattern = re.compile(p1)
        SampleID='#SampleID'
        for ele in self.metadata.columns:
            if pattern.findall(ele):
                if len(pattern.findall(ele)[0]) > 4:
                    SampleID = ele
                    break
            else:
                SampleID='#SampleID'
        self.metadata = self.metadata.set_index(SampleID)
        try:
            self.metadata = self.metadata.drop(['#q2:types'])
        except:
            print('no #q2:types exist')
        #self.merged = metadata.merge(feature_table, left_on=SampleID, right_index=True)
        

    def get_normalize_feature_table(self):
        self.feature_table_bakckup = copy.copy(self.feature_table)
        self.feature_table = self.feature_table.transpose().to_dense()
        for row in self.feature_table.index:
            self.feature_table.loc[row] /= sum(self.feature_table.loc[row])
        self.feature_table = self.feature_table.to_sparse().transpose()
    #def rarefied_feature_table(self,seqs_per_sample):
    #    pass

    def recursion_tree(self,node):
        """recursion to get the sample_series of a tree
        """
        if node.clades: # for non-leaf node
            tmp = 0
            flag = 0
            for clade in node.clades:
                if flag == 0:
                    tmp = copy.copy(self.recursion_tree(clade).sample_series)
                else:
                    tmp += self.recursion_tree(clade).sample_series   
                flag = 1
            node.sample_series = tmp
        else: # leaf node which has been init above.
            try:
                a = node.sample_series
                #print(node.name +' is a leaf')
            except:
                print('please initialize the tree leaves by otu table.')
        return node
    def get_node_domain_otu(self,node):
        """recurse to get the domain otu of every node from the root to leaves."""
        if node.clades:
            if node.clades[0].abu < node.clades[1].abu:
                node.domain_otu = self.get_node_domain_otu(node.clades[1]).domain_otu
                self.get_node_domain_otu(node.clades[0])
            else:
                node.domain_otu = self.get_node_domain_otu(node.clades[0]).domain_otu
                self.get_node_domain_otu(node.clades[1])
        return node
        
                        

    def get_domain_otu(self):
        """ for non  terminal node, we need to get the domain otu."""
        for leaf in self.feature_tree.get_terminals():
            leaf.domain_otu = leaf.name
        self.feature_tree = self.get_node_domain_otu(self.feature_tree)
        #print(temp)

        """
        for clade in self.feature_tree.find_clades(order='level'):
            if not clade.clades:
                continue
            if clade.clades[0].abu < clade.clades[1].abu:
                clade.domain_otu = clade.clades[1].name
            else:
                clade.domain_otu = clade.clades[0].name
        print('self feature tree root abu')
        print(self.feature_tree.domain_otu)
        """


    


    def get_featured_tree(self):
        """ For every node in the tree can be seen as a feature of the data
        (leaves are just OTUs,internal nodes are combinded by OTUs).
        We want to get the feature's data (i.e. columns in OTU table).
        For example,  node i is  the parent of node j and node k which are 
        leaves.So the column of node i and node j can be obtained from the
        OTU table directly.Consequently,The column of node i will be
        generate from the column of node k and node j.

        node_i  = {Sample_0:node_j[Sample_0]+node_k[Sample_0],
                   Sample_1:node_j[Sample_1]+node_k[Sample_1],
                   ...
                   Sample_n:node_j[Sample_n]+node_k[Sample_n],
        }
        """

        for t in self.tree.get_terminals():
            t.sample_series = self.feature_table.transpose()[t.name]
        self.feature_tree = self.recursion_tree(self.tree.root)
        i = 0
        for clade in self.feature_tree.find_clades(order='level'):
            clade.ID_num = i 
            clade.abu = np.mean(clade.sample_series.values)
            #clade.domain_otu = clade.sample_series.idxmax()
            i += 1

    def get_subtree(self,ID_num=0):
        """get the subtree whose root's ID_num equal to the given parameter.
        """
        for clade in self.feature_tree.find_clades(order='level'):
            if ID_num == clade.ID_num:
                self.subtree = clade
                break

    @staticmethod
    def split_features_value(feature_col, label_col):
        """split the feature  value (x1,x2,x3,...,xn) into (x1,..,xi) and
        (x_i+1,xn) according to the label. e.g. Healthy and Disease
        Arg: 
            feature_col: values for samples in a feature.
            label_col: label for samples.
        Return:
            two arrays with different label
        """
        feature = feature_col.sort_index()
        label = label_col.sort_index()
        tmp = ''    
        part1 = []
        part2 = []
        for i in range(len(label)):
            ele = label[i]
            if tmp == '':
                tmp = ele
            if ele == tmp:
                part1.append(feature[i])
            else:
                part2.append(feature[i])
        return part1,part2

    def stats_test(self, label_col, method_name,ID_num):
        """traverse the tree and compute the p value for every node(feature).
        Args:
            col: a Series object  it is key-value pair (SampleID,col_value)
        """
        #split_two_series =
        self.get_subtree(ID_num)
        label_col = self.metadata[label_col]
        label_col = label_col.sort_index()
        methods = {'F_test':stats_test.F_test,
               't_test':stats_test.t_test,
               'fisher_exact_test':stats_test.fisher_exact_test,
               'mannwhitneyu':stats_test.mannwhitneyu_test}
        method = methods[method_name]
        for node in self.subtree.find_clades(order='level'):
            feature = node.sample_series.sort_index()
            part1 ,part2 = self.split_features_value(feature,label_col)
            pvalue = method(part1, part2)
            node.pvalue = -math.log(pvalue)
    
    def get_corr_coefficient(self,label_col, method_name,ID_num):
        self.get_subtree(ID_num)
        label_col = self.metadata[label_col]
        label_col = label_col.sort_index()
        methods = {'spearman':stats.spearmanr,
                   'pearson':stats.pearsonr
        }
        method = methods[method_name]
        for node in self.subtree.find_clades(order='level'):
            feature = node.sample_series.sort_index()
            values = [float(value) for value in feature.values]
            node.corr_coef = method(values,label_col)[0]

    def get_GI(self):
        """Compute Gini Index for every node .
            GI = 1 - \frac{\sum x_i^2}{(\sum x_i)^2}
        """
        for clade in self.feature_tree.find_clades(order='level'):
            total = sum(clade.sample_series.values)**2
            seperate = sum([ele**2 for ele in clade.sample_series.values])
            try:
                clade.GI = 1-(seperate/total)
            except:
                print('error occured when compute Gini Index')
    def plot_scatter(self,para1,para2,ID_num=0):
        """ para1 can be abundance or GI or pvalue or corr_coef
        """
        self.get_subtree(ID_num)
        dict_of_values = {
            'GI':[],
            'pvalue':[],
            'abundance':[],
            'corr_coef':[]
        }
        for clade in self.subtree.find_clades(order='level'):
            try:
                dict_of_values['GI'].append(clade.GI)
            except:
                pass
                #print('no GI in the clade')
            try:
                dict_of_values['abundance'].append(clade.abu)
            except:
                pass
                #print('no abu in the clade')
            try:
                dict_of_values['pvalue'].append(clade.pvalue)
            except:
                pass
                #print('no pvalue in the clade')
            try:
                dict_of_values['corr_coef'].append(clade.corr_coef)
            except:
                pass
                #print('no corr_coef in the clade')
        names = ['ID num:'+ str(clade.ID_num) for clade in self.subtree.\
            find_clades(order='level')]
        colors = [clade.plot_color for clade in self.subtree.find_clades\
            (order='level')]
        traces = []
        for i in range(len(colors)):
            trace = plotly.graph_objs.Scatter(
                x = [dict_of_values[para1][i]],
                y = [dict_of_values[para2][i]],
                name = '',
                mode = 'markers',
                marker = dict(color = colors[i]),
                showlegend = False,
                text = names[i]
            )
            traces.append(trace)
        '''
        trace = plotly.graph_objs.Scatter(
            x = dict_of_values[para1],
            y = dict_of_values[para2],
            mode = 'markers',
            text = names
        )
        '''
        layout = plotly.graph_objs.Layout(
            title = para1 +' and ' +para2 ,
            xaxis =dict(title=para1),
            yaxis =dict(title=para2),
            hovermode = 'closest'
        )
        fig = plotly.graph_objs.Figure(
            data = traces,
            layout = layout
        )
        div = plotly.offline.plot(fig,output_type='div')
        return div
    def plot_tree(self):
        div = rectangle_subtree.plot_tree(self.subtree)
        return div
    def plot_whole_tree(self):
        GI_arr = []
        abu_arr = []
        errror_clade_count = 0
        for clade in self.feature_tree.find_clades(order='level'):
            try:
                #a = clade.GI
                #print('GI exist')
                GI_arr.append(clade.GI)
                abu_arr.append(clade.abu)
            except:
                #print('error when dealing with the whole tree')
                errror_clade_count += 1
        print('eerrr clade count %d '%errror_clade_count)
        names = ['ID num:'+ str(clade.ID_num) for clade in self.feature_tree.\
            find_clades(order='level')]
        colors = [clade.plot_color for clade in self.feature_tree.find_clades\
            (order='level')]
        traces = []
        for i in range(len(colors)):
            trace = plotly.graph_objs.Scatter(
                x = [abu_arr[i]],
                y = [GI_arr[i]],
                name = '',
                mode = 'markers',
                marker = dict(color = colors[i]),
                showlegend = False,
                text = names[i]
            )
            traces.append(trace)
            '''
        trace = plotly.graph_objs.Scatter(
            x = abu_arr,#dict_of_values[para1],
            y = GI_arr,#dict_of_values[para2],
            mode = 'markers',
            text = names
        )
        '''
        layout = plotly.graph_objs.Layout(
            title = ' Abudance and GI ' ,
            xaxis =dict(title='Abudance'),
            yaxis =dict(title='GI'),
            hovermode = 'closest'
        )
        fig = plotly.graph_objs.Figure(
            data = traces,
            layout = layout
        )
        div = plotly.offline.plot(fig,output_type='div')
        return div
    #def colored_tree_plot(self,colors,taxon):
    #    for clade in tree.find_clades()

        






        