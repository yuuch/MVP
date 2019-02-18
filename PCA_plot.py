import numpy as np
import math
import pandas as pd
from plotly import graph_objs as go
import plotly
from sklearn.decomposition import PCA
from plotly import tools
def func1(M,numeric_cols, axises_num=[0,1],n=4):
    """ PCA the matrix, plot the scatter and plot some axises.
        Args: M : numpy.ndarray
            axises_num: the order number of the axises  
    """
    pca = PCA()
    pca.fit(M)
    scatters = pca.fit_transform(M)
    #print(len(scatters))
    pairs = []
    for i in range(min(4,n)):
        for j in range(i):
            pairs.append((j,i))
    traces_array = []
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
        #print(len(y))
    # get axis
        n = len(pca.components_[0])
        axis_results = []
        for i in range(n):
            one_hoc = np.zeros(n)
            one_hoc[i]  = 1
            inv = np.linalg.inv(pca.components_)
            single_result = np.dot(inv,one_hoc)
            axis_results.append(single_result)
        traces = []
        trace0 = go.Scatter(x=x,y=y,name='',mode='markers')
        traces.append(trace0)
    # adjust the length of project axises
        for i in axises_num:
            temp0 = 0
            temp1 = 0
            temp3 = axis_results[i][first]**2+axis_results[i][second]**2 
            if temp3 < Max:
                temp0 = math.sqrt(Max/temp3)*axis_results[i][first]
                temp1 = math.sqrt(Max/temp3)*axis_results[i][second]
            trace = go.Scatter(x=[0,temp0],y = [0,temp1],mode='lines',name = numeric_cols[i])
            traces.append(trace)
        traces_array.append(traces)
    return traces_array,pairs
    
def six_subplot(traces_array, pairs, num=6):
    fig = tools.make_subplots(rows=2,cols=3)
    for i in range(len(traces_array)):
        if i < 3:
            for trace in traces_array[i]:
                fig.append_trace(trace,1,i+1)
        else:
            for trace in traces_array[i]:
                fig.append_trace(trace,2,i+1-3)
    for i in range(len(pairs)):
        fig['layout']['xaxis'+str(i+1)].update(title='PC'+str(pairs[i][0]))
        fig['layout']['yaxis'+str(i+1)].update(title='PC'+str(pairs[i][1]))
    div = plotly.offline.plot(fig,output_type='div')
    return div




    pass 
def judge_numeric_col(df):
    num_cols = []
    for col in df.columns:
        try:
            int(df[col].values[0])
            num_cols.append(col)
        except:
            pass
    return num_cols
def run_this_script(df):
    numeric_cols = judge_numeric_col(df)
    #print(numeric_cols)
    M = df[numeric_cols].values
    traces_array ,pairs = func1(M, numeric_cols)
    div = six_subplot(traces_array, pairs)
    return div
if __name__ == "__main__":
    df = pd.read_csv('upload_files/demo_metadata.tsv', sep='\t')
    #### drop the row named #q2 type
    df = df.drop(0)

    div = run_this_script(df)
    f = open('Feb15.html', 'w')
    f.write(div)