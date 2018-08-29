from Bio import Phylo
from plotly.offline import download_plotlyjs, init_notebook_mode,  iplot, plot
init_notebook_mode(connected=True)
import plotly.graph_objs as go
def get_y_coordinates(tree):
    """
    Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    #print(xcoords)
    #tree.depth() maps tree clades to depths (by branch length).
    #returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth is the distance from root
    #to clade
    
    #  If there are no branch lengths, assign unit branch lengths
    
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    else:
        M = max(xcoords.values())
        for ele in xcoords:
            xcoords[ele]=M-xcoords[ele]
    #print('xcoords')
    #print(xcoords)
    return xcoords
   
def get_x_coordinates(tree, dist=1):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
        """
    maxheight = tree.count_terminals()#Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i*dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
    def calc_row(clade):
            for subclade in clade:
                if subclade not in ycoords:
                    calc_row(subclade)
            ycoords[clade] = (ycoords[clade.clades[0]] +
                              ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords
#traces=[]
def get_lines(tree,clade,xcoords,ycoords,traces):
    path = tree.get_path(clade)
    x0 = xcoords[clade]
    y0 = ycoords[clade]
    if len(path)>1:
        parent= path[-2]
        x_p = xcoords[parent]
        y_p = ycoords[parent]
        trace_p = go.Scatter(# trace parent vertical line
            x = (x0,x0),
            y = (y0,y_p),
            marker=dict(color='rgb(25,25,25)'),
            mode = 'lines',
            hoverinfo=None
        )
        traces.append(trace_p)
    else:
        if len(path) == 1: # this calde is the subclade of root
            y_r = ycoords[tree.root]
            trace_p = go.Scatter(x=(x0,x0),
                                y=(y0,y_r),
                                 mode='lines',
                                 marker =dict(color='rgb(25,25,25)'),
                                 hoverinfo=None
                                )
            traces.append(trace_p)
            
    if clade.clades:
        x_l = xcoords[clade.clades[0]]
        x_r = xcoords[clade.clades[-1]]
        trace_c = go.Scatter(
        x = (x_l,x_r),
        y = (y0,y0) ,
        marker=dict(color='rgb(25,25,25)'),
        mode = 'lines',
        hoverinfo=None
        )
        traces.append(trace_c)
        for ele in clade.clades:
            get_lines(tree,ele,xcoords,ycoords,traces=traces)
    else:
        if y0 >0:
            trace_t = go.Scatter(
            line = dict(color='rgb(25,25,25)',dash='dash'),
            mode = 'lines',
            x =(x0,x0),
            y=(y0,0))
            traces.append(trace_t)
def plot_tree(tree):
    text_dict ={}
    for clade in tree.find_clades(order='level'):
        if clade.name:
            text_dict[clade]='<br>name:'+clade.name+'<br> clade num:'+str(clade.num)
        else:
            text_dict[clade]='<br>clade num:'+str(clade.num)
    #tree = Phylo.read('test_tree.nwk','newick')
    y_coords = get_y_coordinates(tree)
    x_coords = get_x_coordinates(tree)
    X = []
    Y = []
    text = []
    for key in x_coords.keys():
        X.append(x_coords[key])
        Y.append(y_coords[key])
        text.append(text_dict[key])
        #words = key.name
        '''
        if hasattr(key,'name'):
            words = 'name: '+key.name+'<br>node_num'+str(key.node_num)
        else:
            words = '<br>node_num'+str(key.node_num)
        '''
        #text.append(words)
    traces = []
    get_lines(tree,tree.root,x_coords,y_coords,traces=traces)
    print('len straces')
    print(len(traces))
    trace = go.Scatter(
        x = X,
        y = Y,
        mode = 'markers',
        text= text,
        name=''
    )
    #print(X)
    #print(Y)
    data = [trace]
    layout = dict(showlegend=False)
    for ele in traces:
        data.append(ele)
    fig =go.Figure(data=data,layout=layout)
    tree_div = plot(fig,output_type='div')
    return tree_div
