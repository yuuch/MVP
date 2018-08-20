import biom
import pandas as pd
import copy
import numpy as np
import plotly
class Heatmap(object):
    def __init__(self,meta_data, feature_table, file_type='biom'): #metadata is sample metadata
        self.meta = pd.read_csv(meta_data,sep='\t')
        self.biom_table = biom.load_table(feature_table)
        self.df = self.biom_table.to_dataframe().transpose()
        self.df_primary_col = self.df.columns

    def map(self):
        """
        map feature-table to sample-metadata, merge two table
        """
        self.df_primary_col = self.df.columns
        array1 = self.meta[self.meta.columns[0]]# sample id of metadata
        array2 = self.df[self.df.columns[0]]  # sample id of feature-table
        mapped_dict ={'metadata':[],'feature_table':[]}
        for i in range(len(array1)):
            for j in range(len(array2)):
                if array2[j] == array1[i]:
                    mapped_dict['metadata'].append(i)
                    mapped_dict['feature_table'].append(j)
                    break
        print(mapped_dict['feature_table'])
        temp_table = self.df.iloc[mapped_dict['feature_table'],:]
        temp_table.index = list(range(temp_table.shape[0]))
        temp_meta = self.df.iloc[mapped_dict['metadata'],:]
        temp_meta.index = list(range(temp_meta.shape[0]))
        assert temp_meta.shape[0] == temp_table.shape[0]
        self.df = pd.concat([temp_table,temp_meta],axis=1)

    def sort_by_features(self,*args):
        list_args = [ele for ele in args]
        if ''in list_args:
            list_args.remove('')
        if 'None' in list_args:
            list_args.remove('None')
        if len(list_args)>0:
            self.df = self.df.sort_values(by=list_args)

    def obtain_numerical_matrix(self):
        new_df = self.df[self.df_primary_col]
        self.df = new_df

    def filter(self,prevalence_threshold=0.1, abundance_num=100, variance_num=100):
        """
        ignor low prevalence ,low abundance feature and low variance.
        choose the first 100 abundance feature and varivance feature.
        """
        temp_col = copy.copy(self.df_primary_col)
        temp_df = self.df[temp_col]

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
    def plotly_div(self):
        '''
        the next three lines are test df.
        self.df = pd.DataFrame({
            'col0':[0,1,2,3,4,5],
            'col1':[1,2,3,4,5,6]
        })
        '''
        data = [plotly.graph_objs.Heatmap(z=self.df.values.tolist())]
        layout = plotly.graph_objs.Layout(
            width = 1000,
            height = 800
        )
        fig = plotly.graph_objs.Figure(data=data,layout=layout)
        div_string = plotly.offline.plot(fig,filename='plotly.html',output_type='div')
        return div_string


    