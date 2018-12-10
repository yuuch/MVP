import biom
import os
import numpy as np
import pandas as pd
from skbio.stats import subsample_counts
import skbio
import plotly 
from plotly import graph_objs as go
import re
def get_rarefied(otu_table,seqs_per_sample):
    """
    Args:
        otu_table:(dataframe) load biom file and change to dataframe
        seqs_per_sample:...
    Rerutn:
        a rarefied OTU table
    """
    new_counts = []
    for sample in otu_table.columns:
        arr = []
        seqs = sum(otu_table[sample])
        if seqs <= seqs_per_sample:
            arr = np.array(otu_table[sample].values).astype(int)
        else:
            values = np.array(otu_table[sample].values).astype(int)
            arr = subsample_counts(values,seqs_per_sample)
        new_counts.append(arr)
    rarefied = pd.DataFrame(new_counts,columns=otu_table.index,index=otu_table.columns)
    return rarefied.T

def alpha_diversity_compute(df,metric,tree=None):
    #df = biom.load_table(otu_table).to_dataframe()
    result = ''
    if metric == 'faith_pd':
        tree = skbio.TreeNode.read(tree)
        result = skbio.diversity.alpha_diversity(
            counts=df.T.values,ids=df.columns,
            metric='faith_pd',tree=tree,
            otu_ids=df.index)
    elif metric == 'ace':
        result = skbio.diversity.alpha_diversity(counts=df.T.values.astype(int), ids=df.columns,metric=metric)
    else:
        result = skbio.diversity.alpha_diversity(counts=df.T.values, ids=df.columns,metric=metric)
    #result = pd.DataFrame(result,columns=['alpha_div'])
    return result

def get_rarefied_names(rarefy_result,rarefied_num=10):
    names = rarefy_result.columns
    p1 = '-.*?-'
    pattern = re.compile(p1)
    result_dict  = {}
    for name in names:
        seqs_num = pattern.findall(name)[0][1:-1]
        if seqs_num in result_dict:
            result_dict[seqs_num].append(name)
        else:
            result_dict[seqs_num]=[name]
    return result_dict
#grouped_columns = get_rarefied_names(rarefy_result)   

def alpha_box_plot(box_data):
    x = []
    y = {}
    for feature in box_data:
        tmp = box_data[feature]
        y[feature]=[]
        for depth in tmp:
            for ele in tmp[depth]:
                x.append(depth)
                y[feature].append(ele)
    data =[]
    for feature in y:
        trace = go.Box(
            y=y[feature],
            x=x,
            name=feature)
        data.append(trace)
    layout = go.Layout(boxmode='group')
    fig = go.Figure(data=data,layout=layout)
    return plotly.offline.plot(fig,output_type='div')
#alpha_box_plot(result)

def get_rarefy_result(otu_table_path, seq_max, step, metric, tree=None,rarefied_num=10):
    seqs_array = list(range(step,seq_max,step))
    #print(seqs_array)
    otu_table = biom.load_table(otu_table_path).to_dataframe()
    #seqs_alpha_dict = {}
    result = pd.DataFrame([],index=otu_table.columns)
    for seqs_per_sample in seqs_array:
        #print(result.columns)
        for i in range(rarefied_num):
            rarefied = get_rarefied(otu_table,seqs_per_sample)
            rarefied_alpha =alpha_diversity_compute(rarefied,metric,tree=tree)
            rarefied_df = pd.DataFrame(rarefied_alpha,columns=['depth-'+str(seqs_per_sample)+'-iter-'+str(i)])
            result = result.merge(rarefied_df,left_index=True,right_index=True)
    return result

def plot_alpha_rarefaction(feature_table,metadata,max_seq,step,metric,meta_column,tree=None):
    rarefy_result = get_rarefy_result(feature_table,max_seq,step,metric,tree=tree)
    metadata = pd.read_csv(metadata,sep='\t')
    p1 = '.*[Ss][Aa][Mm][Pp][Ll][Ee].*[Ii][Dd].*'
    pattern = re.compile(p1)
    SampleID='#SampleID'
    for ele in metadata.columns:
        if len(pattern.findall(ele)[0]) > 5:
            SampleID = ele
            break
        else:
            SampleID='#SampleID'
    metadata.set_index(metadata[SampleID])
    merged = rarefy_result.merge(metadata,left_index=True, right_on=SampleID)
    grouped_columns = get_rarefied_names(rarefy_result) 
    result = get_box_plot_data(merged,grouped_columns,meta_column=meta_column)
    return alpha_box_plot(result)
#str1 = plot_alpha_rarefaction('feature-table.biom','demo_metadata.tsv',500,40,'observed_otus','BodySite',tree='tree.nwk')
#f = open('alpha_rarefaction.html','w')
#f.write(str1)

def get_rarefied_names(rarefy_result,rarefied_num=10):
    names = rarefy_result.columns
    p1 = '-.*?-'
    pattern = re.compile(p1)
    result_dict  = {}
    for name in names:
        seqs_num = pattern.findall(name)[0][1:-1]
        if seqs_num in result_dict:
            result_dict[seqs_num].append(name)
        else:
            result_dict[seqs_num]=[name]
    return result_dict
#grouped_columns = get_rarefied_names(rarefy_result)   

#merged =merged.set_index('#SampleID')

def get_box_plot_data(merged_df, grouped_columns ,meta_column='BodySite'):
    """get data for box plot.
    Args:
        merged_df: mrege otutable with metadata.
        grouped_columns: a dict keep columns with the same seqs-depth.
            e.g   {'400':['depth-400-iter-0',...,'depth-400-iter-9']}
    Return:
        a dict contain meta column and seqs-depth information.
        e.g. {'gut':{'400':[1,2,3,4]},
              'left_palm':{'400':[1,3,4,5]}
              }
    """
    #old = pd.DataFrame([])
    appeared_features = []
    dfs = {}
    for ele in merged_df[meta_column]:
        if ele in appeared_features:
            continue
        else:
            appeared_features.append(ele)
        mask = merged_df[meta_column] == ele
        dfs[ele] = merged_df[mask]
    result = {}
    for sample in dfs:
        result[sample]={}
        for col in grouped_columns:
            tmp = dfs[sample][grouped_columns[col]]
            tmp = tmp.values.flatten()
           #print(tmp.shape)
            result[sample][col]=tmp
            #print(tmp)   
        #break
    return result
#result = get_box_plot_data(merged,grouped_columns)
#result['left palm']['400'].shape