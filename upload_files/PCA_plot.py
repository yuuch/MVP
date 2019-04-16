import numpy as np
import math
from plotly import graph_objs as go
import plotly
from sklearn.decomposition import PCA
def func1(M,axises_num=[0,1,2]):
    """ PCA the matrix, plot the scatter and plot some axises.
        Args: M : numpy.ndarray
            axises_num: the order number of the axises  
    """
    pca = PCA()
    pca.fit(M)
    scatters = pca.fit_transform(M)
    # get scatters
    x = []
    y = []
    Max = 0
    for ele in scatters:
        temp = ele[0]**2 + ele[1]**2
        if Max < temp:
            Max = temp
        x.append(ele[0])
        y.append(ele[1])
    # get axis
    n = len(pca.components_[0])
    axis_results = []
    for i in range(n):
        zeros = np.zeros(n)
        zeros[i]  = 1
        inv = np.linalg.inv(pca.components_)
        single_result = np.dot(inv,zeros)
        axis_results.append(single_result)
    # adjust the length of project axises
    traces = []
    trace0 = go.Scatter(x=x,y=y,name='',text='scatters',mode='markers')
    traces.append(trace0)
    for i in axises_num:
        temp0 = 0
        temp1 = 0
        temp3 = axis_results[i][0]**2+axis_results[i][1]**2 
        if temp3 < Max:
            temp0 = math.sqrt(Max/temp3)*axis_results[i][0]
            temp1 = math.sqrt(Max/temp3)*axis_results[i][1]
            
        trace = go.Scatter(x=[0,temp0],y = [0,temp1],mode='lines')
        traces.append(trace)
    div = plotly.offline.plot(traces,output_type='div')
    return div
