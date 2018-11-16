import plotly
import plotly.graph_objs as go
import sklearn.manifold
def dimension_reduce_visualize(X,labels,flag_3d=False):
    """ visualize the result of dimension reduction by ISOmap or MDS
    Args:
         X: the coordinates of the scatter.
         labels: decide the color of the points.
         flag_3d: False for 2D graph ,True for 3d graph.
    Return: a div format string.
    """
    labels_dict = {}
    for i in range(len(labels)):
        if labels[i] in labels_dict:
            labels_dict[labels[i]].append(i)
        else:
            labels_dict[labels[i]] = [i]
    traces = []
    for label in labels_dict:
        if flag_3d == 'False':
            trace = go.Scatter(
                x = X[labels_dict[label], 0],
                y = X[labels_dict[label], 1],
                mode = 'markers',
                name = str(label)
            )
            traces.append(trace)
        else:
            trace = go.Scatter3d(
                x = X[labels_dict[label], 0],
                y = X[labels_dict[label], 1],
                z = X[labels_dict[label], 2],
                mode = 'markers',
                name = str(label)
            )
            traces.append(trace)
    layout = go.Layout(
        title = 'dimension reduction'
    )
    fig = go.Figure(data=traces, layout =layout)
    div_str = plotly.offline.plot(traces, output_type='div')
    return div_str
def reduce_dimension(matrix, n_component, method='Isomap'):
    """ reduce a n*m matrix to n*2 or n*3 matrix.
    Args:
        method: Isomap or MDS.
        n_component: usually 2 or 3.
        matrix: n*m origin matrix
    Return:
        return a n*2 or n*3 matrix
    """
    methods = {
        'Isomap':sklearn.manifold.Isomap,
        'MDS':sklearn.manifold.MDS
    }
    method_class = methods[method]
    embedding = method_class(n_component)
    result = embedding.fit_transform(matrix)
    return result


    
        
