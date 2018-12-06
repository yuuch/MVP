import pandas as pd
import plotly
from plotly import graph_objs as go
import sklearn.manifold
import numpy as np
import skbio
import biom
import scipy
def alpha_diversity_pre(otu_table,metric,tree=None):
    df = biom.load_table(otu_table).to_dataframe()
    result = ''
    if metric == 'faith_pd':
        tree = skbio.TreeNode.read(tree)
        result = skbio.diversity.alpha_diversity(
            counts=df.T.values,ids=df.columns,
            metric='faith_pd',tree=tree,
            otu_ids=df.index)
    else:
        result = skbio.diversity.alpha_diversity(counts=df.T.values, ids=df.columns,metric=metric)
    result = pd.DataFrame(result,columns=['alpha_div'])
    return result
def alpha_diversity(alpha_table, metadata,label_col):
    """ split the alpha table into serveral parts according to the metadata.
    Args:
        alpha_table: an pandas dataframe which come from the 'alpha_dive
            rsity_pre' function.
        metadata: record the macro feature of every sample.(File name, String)
    Return :
        a dict contains every label and its samples.
        e.g. dict0 = {'class0':[0,1,2,3,4,5], 'class1': [5, 6, 7, 8, 9]}
    """
    metadata = pd.read_csv(metadata, sep='\t')
    #alpha_table = pd.read_csv(alpha_table, sep='\t')
    try:
        merged = alpha_table.merge(
            metadata, left_index=True, right_on='#SampleID')
    except:
        print('Wrong column name')
    diversity = merged['alpha_div']
    labels = merged[label_col]
    result_dict = {}
    for j in range(len(labels)):
        i = j+1
        key = labels[i]
        if key  in result_dict:
            result_dict[key].append(diversity[i])
        else:
            result_dict[key] =[diversity[i]]
    return result_dict

def alpha_box_plot(result_dict):
    data = []
    for ele in result_dict:
        tmp_str = '(n='+str(len(ele))+')'
        trace = go.Box(
            name = ele+tmp_str,
            y = result_dict[ele]
        )
        data.append(trace)
    layout = go.Layout(
        title = "alpha diversity"
    )
    fig = go.Figure(data=data, layout=layout)
    div = plotly.offline.plot(fig,output_type='div')
    return div

def Bray_Curtis_distance(otu_table):
    """compute pairwise distances of samples.
    Args:
        otu_table: otu_table file name (biom format)
    Return:
        pairwise distances in DataFrame format
    """
    otu_table = biom.load_table(otu_table).to_dataframe()
    samples = otu_table.columns
    M = np.zeros(shape=(len(samples),len(samples)))
    for i in range(len(samples)):
        for j in range(len(samples)):
            if i == j:
                M[i][j]=0
            elif i>j:
                M[i][j] = M[j][i]
            else:
                M[i][j]= scipy.spatial.distance.braycurtis(otu_table[samples[i]],otu_table[samples[j]])
    df = pd.DataFrame(M,columns=samples,index=samples)
    return df  
    
def beta_diversity_pre(otu_table, tree=None, metric='weighted_unifrac'):
    try: # beta divesity related to the phylo tree
        tree = skbio.TreeNode.read(tree)
        df = biom.load_table(otu_table).to_dataframe()
        unifrac = skbio.diversity.beta_diversity(
            counts=df.T.values, ids=df.columns,metric=metric,
            tree=tree,otu_ids=df.index)
        distance_matrix = pd.DataFrame(
            unifrac.data,
            columns=unifrac.ids,
            index=unifrac.ids
            )
    except:  # do not need the tree
        distance_matrix = Bray_Curtis_distance(otu_table)
    return distance_matrix

def beta_diversity(distance_matrix, metadata_file, n_components=2, col='BodySite',dim_method='Isomap'):
    """ obtain the visualize of the distance matrix.
    Args:
        distance_matrixea:(dataframe)
            distance between samples come frome the beta_diversity_pre
            function
    Return:
        a dict storing the points from the same label.
    """
    metadata = pd.read_csv(metadata_file,sep='\t')
    #df = pd.read_csv(distance_matrix, sep='\t')
    #df = df.set_index(df.columns[0])
    df = distance_matrix
    values= df.values
    # TODO edit Isomap or MDS etc.
    methods = {
        'PCoA': skbio.stats.ordination.pcoa,
        'Isomap': sklearn.manifold.Isomap,
        'MDS': sklearn.manifold.MDS
    }
    method = methods[dim_method]
    try: # manifold method
        embedding = method(n_components=n_components)
        X = embedding.fit_transform(values)
        cols = ['x0','x1']
        if n_components == 3:
            cols.append('x2')
    except: # pcoa method
        dm =skbio.stats.distance.DistanceMatrix(values,ids=df.columns)
        pcoa_result = method(dm,'fsvd',n_components)
        X = pcoa_result.samples.values
        cols = pcoa_result.samples.columns
    value_df = pd.DataFrame(X,index=df.index, columns=cols)
    merged = value_df.merge(metadata,left_index=True,right_on='#SampleID')
    labels = merged[col]
    indexes = merged.index
    coordinates = merged[cols]
    result_dict = {}
    for i in range(len(labels)):
        key = labels.iloc[i]
        if key in result_dict:
            result_dict[key].append(coordinates.iloc[i])
        else:
            result_dict[key] = [coordinates.iloc[i]]
    return result_dict

def plot_beta_scatter(result_dict):
    data = []
    for ele in result_dict:
        tmp = np.array(result_dict[ele])
        if len(result_dict[ele][0]) == 2:
            trace = go.Scatter(
                x = tmp[:,0],
                y = tmp[:,1],
                name = ele,
                mode = 'markers'
            )
            data.append(trace)
        else:
            trace = go.Scatter3d(
                x = tmp[:,0],
                y = tmp[:,1],
                z = tmp[:,2],
                mode = 'markers',
                name = ele
            )
            data.append(trace)
    layout = go.Layout(title="beta diversity")
    fig = go.Figure(data=data, layout=layout)
    div = plotly.offline.plot(fig,output_type='div')
    return div



    







        




    

    
    


