import numpy as np
import math
import pandas as pd
from plotly import graph_objs as go
import plotly
from sklearn.decomposition import PCA
from plotly import tools
import corr_tree_new
from sklearn import preprocessing
def judge_numeric_col(df):
    num_cols = []
    for col in df.columns:
        try:
            int(df[col].values[2])
            num_cols.append(col)
        except:
            pass
    return num_cols


def pca_function(mvp_tree, ID_num,axis_names=None,n=4):

    """ PCA the matrix, plot the scatter and plot some axises.
        Args: M : numpy.ndarray
            axises_num: the order number of the axises  
    """
    numeric_cols = judge_numeric_col(mvp_tree.metadata)
    M = mvp_tree.metadata[numeric_cols]
    M = preprocessing.scale(M)
    pca = PCA()
    pca.fit(M)
    scatters = pca.fit_transform(M)
    #print(len(scatters))
    ratio = pca.explained_variance_ratio_
    pairs = []
    for i in range(min(4,n)):
        for j in range(i):
            pairs.append((j,i))
    traces_array = []
    sample_names = mvp_tree.metadata.index.values
    sample_abuns = mvp_tree.subtree.sample_series
    for pair in pairs:
        first = pair[0]
        second = pair[1]
    # get scatters
        x = []
        y = []
        Max = 0
        for ele in scatters:
            temp = ele[first]**2 + ele[second]**2
            if Max < temp:
                Max = temp
            x.append(ele[first])
            y.append(ele[second])
    # get axis
        n = len(pca.components_[0])
        axis_results = []

        # get the abundance for the clade in the subtree.(mean abundance for all sample)
        colors = [sample_abuns[name] for name in sample_names]
        text = ['abun: '+str(sample_abuns[name]) for name in sample_names]
        for i in range(n):
            one_hoc = np.zeros(n)
            one_hoc[i]  = 1
            inv = np.linalg.inv(pca.components_)
            single_result = np.dot(inv,one_hoc)
            axis_results.append(single_result)
        traces = []
        trace0 = go.Scatter(x=x,y=y,name='',mode='markers',marker=dict(color=colors,\
            colorbar=dict(title='abundance'),colorscale='Viridis'),showlegend=False,text=text)
        traces.append(trace0)

    # adjust the length of project axises
        if axis_names == None:
            axis_names = numeric_cols[0:2]
    # get name order in the numeric cols
        name_order = [numeric_cols.index(ele) for ele in axis_names]
        for i in name_order:
            temp0 = 0
            temp1 = 0
            temp3 = axis_results[i][first]**2+axis_results[i][second]**2 
            if temp3 < Max:
                temp0 = math.sqrt(Max/temp3)*axis_results[i][first]
                temp1 = math.sqrt(Max/temp3)*axis_results[i][second]
            trace = go.Scatter(x=[0,temp0],y = [0,temp1],mode='lines',line=dict(dash='dash'),name = numeric_cols[i],\
                showlegend=False)
            traces.append(trace)
        traces_array.append(traces)
    return traces_array, pairs, ratio
    
def six_subplot(traces_array, pairs, ratio, num=6):
    fig = tools.make_subplots(rows=2,cols=3)
    for i in range(len(traces_array)):
        if i < 3:
            for trace in traces_array[i]:
                fig.append_trace(trace,1,i+1)
        else:
            for trace in traces_array[i]:
                fig.append_trace(trace,2,i+1-3)
    for i in range(len(pairs)):
        fig['layout']['xaxis'+str(i+1)].update(title='PC'+str(pairs[i][0])+': '\
            +"{0:.2f}".format(ratio[pairs[i][0]]*100) + '%')
        fig['layout']['yaxis'+str(i+1)].update(title='PC'+str(pairs[i][1])+': '\
            +"{0:.2f}".format(ratio[pairs[i][1]]*100) + '%') 
    div = plotly.offline.plot(fig,output_type='div')
    return div

def run_this_script(mvp_tree, ID_num,axis_names=None):
    #  df is metadata df ,obtain the numerical cols
    #print(numeric_cols)
    traces_array ,pairs, ratio= pca_function(mvp_tree, ID_num,axis_names)
    div = six_subplot(traces_array, pairs,ratio)
    return div

if __name__ == "__main__":
    df = pd.read_csv('upload_files/demo_metadata.tsv', sep='\t')
    #### drop the row named #q2 type
    df = df.drop(0)
    feature_table_path='upload_files/feature-table.biom'
    tree_path='upload_files/tree.nwk'
    metadata_path='upload_files/demo_metadata.tsv'
    ID_num=114
    mvp_tree = corr_tree_new.MvpTree(feature_table_path,tree_path,metadata_path,ID_num=ID_num)
    div = run_this_script(mvp_tree, ID_num,axis_names=['Year','Day'])
    f = open('Feb15.html', 'w')
    f.write(div)