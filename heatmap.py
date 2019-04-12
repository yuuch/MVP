import biom
import pandas as pd
import copy
import numpy as np
import plotly
import re
from Bio import Phylo

class Heatmap(object):
    def __init__(self,meta_data, feature_table, file_type='biom'): #metadata is sample metadata
        self.meta = pd.read_csv(meta_data,sep='\t')
        self.select_features = self.meta.columns[:3]
        self.biom_table = biom.load_table(feature_table)
        self.df = self.biom_table.to_dataframe().transpose().to_dense()
        self.df = self.df.div(self.df.sum(axis=1), axis=0)
        self.df_primary_col = self.df.columns
        self.numeric_df = copy.copy(self.df)
        self.metadata_set_index()
        #self.ordered = False
        '''
        Here df is a table like:
                otu0 otu1 otu2
        sample0
        sample1
        sample2
        '''
        #print(self.df.index)]
    

    def metadata_set_index(self):
        #metadata = self.metadata
        #feature_table = self.feature_table.transpose()
        p1 = '.*[Ss][Aa][Mm][Pp][Ll][Ee].*[IiNn][DdAa]'
        pattern = re.compile(p1)
        SampleID='#SampleID'
        for ele in self.meta.columns:
            if pattern.findall(ele):
                if len(pattern.findall(ele)[0]) > 4:
                    SampleID = ele
                    break
            else:
                SampleID='#SampleID'
        self.meta = self.meta.set_index(SampleID)
        try:
            self.meta = self.meta.drop(['#q2:types'])
        except:
            print('no #q2:types exist')

    def map(self):
        """
        map feature-table to sample-metadata, merge two table
        """
        self.df_primary_col = self.df.columns
        #print(len(self.df.columns))
        array1 = self.meta.index# sample id of metadata
        array2 = self.df.index  # sample id of feature-table
        mapped_dict ={'metadata':[],'feature_table':[]}
        for i in range(len(array1)):
            for j in range(len(array2)):
                if array2[j] == array1[i]:
                    mapped_dict['metadata'].append(i)
                    mapped_dict['feature_table'].append(j)
                    break

        temp_table = self.df.iloc[mapped_dict['feature_table'],:]
        temp_table.index = list(range(temp_table.shape[0]))
        temp_meta = self.meta.iloc[mapped_dict['metadata'],:]
        temp_meta.index = list(range(temp_meta.shape[0]))
        assert temp_meta.shape[0] == temp_table.shape[0]
        self.df = pd.concat([temp_table,temp_meta],axis=1)
        new_index = []
        for ele in mapped_dict['metadata']:
            new_index.append(array1[ele])
        self.df.index=new_index

    def sort_by_features(self,*args):

        list_args = [ele for ele in args]
        self.select_features = list_args
        copy_args = copy.copy(list_args)
        for ele in copy_args:
            if ele =='':
                list_args.remove('')
            if ele == 'None':
                list_args.remove('None')
        if len(list_args)>0:
            print(list_args)
            self.df = self.df.sort_values(by=list_args)
            """
            self.ordered = True
            self.selected_metadata_df = self.df[list_args]
            """

    def obtain_numerical_matrix(self,cols=None):
        #cols = [ele.name for ele in tree.get_terminals()]
        if not cols:
            cols = self.df_primary_col
            print('cols num: ',len(cols))
        new_df = self.df[cols]
        self.numeric_df = new_df

    def filter(self,prevalence_threshold=0.1, abundance_num=100, variance_num=100):
        """
        ignor low prevalence ,low abundance feature and low variance.
        choose the first 100 abundance feature and varivance feature.
        """
        temp_col = list(copy.copy(self.df_primary_col))
        temp_df = self.df[temp_col]
        prevalence_threshold = float(prevalence_threshold)
        abundance_num = int(abundance_num)
        variance_num = int(variance_num)
        #prevalence 
        for col in self.df_primary_col:
            non_zero_count= 0
            for ele in temp_df[col]:
                if ele > 0:
                    non_zero_count += 1
            prevalence_of_col = non_zero_count/len(temp_df[col])
            if prevalence_of_col < prevalence_threshold:
                temp_col.remove(col)
        # abundance
        abundance_dict = {}
        for col in temp_col:
            abundance = 0
            for ele in temp_df[col]:
                if ele > 0:
                    abundance += ele
            abundance_dict[col] =  abundance

        if len(abundance_dict) < abundance_num:
            pass
        else:
             # sort and select abundance_num features
             sort_abun = sorted(abundance_dict,key=abundance_dict.get,reverse=True)
             temp_col = sort_abun[0:abundance_num]
        # variance
        variance_dict = {}
        for col in temp_col:
            std_of_col = np.std(temp_df[col])
            variance_dict[col] = std_of_col
        if len(variance_dict) < variance_num:
            pass
        else:
            sort_vari = sorted(variance_dict,key=variance_dict.get,reverse=True)
            temp_col = sort_vari[0:variance_num]
        self.df = self.df[temp_col]

    @staticmethod
    def map_metadata_values(selected_metadata_df):
        """map a matadata df(string value) to a number value.
        and save the origin values in the hovertext obj.
        """
        metadata_df = copy.copy(selected_metadata_df)
        hovertext =selected_metadata_df.values
        appeared = {}
        i = 1
        for col in metadata_df.columns:
            for j,ele in enumerate(metadata_df[col]):
                if ele in appeared:
                    metadata_df[col][j] =appeared[ele]
                else:
                    appeared[ele] = i
                    metadata_df[col][j] = i
                    i += 1
        return metadata_df, hovertext

    def plotly_div(self,show_metadata_label=True):
        '''
        the next three lines are test df.
        self.df = pd.DataFrame({
            'col0':[0,1,2,3,4,5],
            'col1':[1,2,3,4,5,6]
        })
        '''
        trace1 = plotly.graph_objs.Heatmap(z=self.numeric_df.values.tolist(),
                                        x=list(self.df.columns),y=self.df.index,
                                        xaxis='x2',colorscale='Viridis')
        trace2 = plotly.graph_objs.Heatmap()
        select_metadata = self.df[self.select_features]
        
        metadata_df, hovertext = self.map_metadata_values(select_metadata)
        trace2 = plotly.graph_objs.Heatmap(z=metadata_df.values,
                 x = list(metadata_df.columns), y = self.df.index,
                 text=hovertext,hoverinfo='text')
        data = [trace2,trace1]
        try:
            assert show_metadata_label
            layout = plotly.graph_objs.Layout(xaxis=dict(
                    domain=[0,0.1]
                ),xaxis2=dict(
                    domain = [0.10001,1]
                )
                )   
        except:
            layout = plotly.graph_objs.Layout()
        fig = plotly.graph_objs.Figure(data=data,layout=layout)
        div_string = plotly.offline.plot(fig,filename='plotly.html',output_type='div')
        return div_string


    