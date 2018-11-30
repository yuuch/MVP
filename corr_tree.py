from Bio import Phylo
from plotly.offline import download_plotlyjs, init_notebook_mode,  iplot, plot
init_notebook_mode(connected=True)
import plotly.graph_objs as go
import math
import gini_index_compute
import scipy.stats
import pandas as pd
import stats_test
import biom
import numpy as np
def get_branch_depth(clade):
    if clade.clades:
        clade.child ={}
        for subclade in clade.clades:
            #print(subclade.name)
            clade.child[subclade]= get_branch_depth(subclade)+1
        clade.depth = max(clade.child.values())
        return max(clade.child.values())
    else:
        clade.child = {}
        clade.depth = 0
        return 0
def get_y_coordinates(tree):
    """
    Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    '''
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
    #M = max(xcoords.values())
    #for ele in xcoords:
    #    xcoords[ele]=math.log2(1+((xcoords[ele]+M)/(M+M)))
    #    xcoords[ele]=-math.exp(-(xcoords[ele])**2)
    '''
    """
    xcoords = tree.depths(unit_branch_lengths=True)
    M = max(xcoords.values())
    for ele in xcoords:
        xcoords[ele]=M-xcoords[ele]
    """
    '''
    for ele in tree.get_terminals():
        xcoords[ele]=0
    '''
    xcoords = {}
    get_branch_depth(tree.root)
    for clade in tree.find_clades(order='level'):
        xcoords[clade] = clade.depth
    return xcoords

def get_phylum_colors(taxonomy_file):
    colors  = ['rgb(255,0,0)','rgb(255,247,0)','rgb(255,0,247)','rgb(162,255,0)',
        'rgb(255,111,0)','rgb(111,0,255)','rgb(2,104,23)','rgb(104,2,40)','rgb(2,23,104)',
        'rgb(192,77,11)','rgb(12,192,162)','rgb(126,3,145)','rgb(17,102,7)','rgb(57,12,179)',
        'rgb(195,216,8)','rgb(216,8,154)','rgb(6,97,94)','rgb(97,6,36)','rgb(125,218,5)']
    df = pd.read_csv(taxonomy_file,sep='\t')
    dict1 = {}
    for i in range(len(df['Taxon'])):
        tmp_str = df.loc[i]['Taxon'].split(';')
        try:
            key = tmp_str[1]
        except:
            key = tmp_str[0]
        if key in dict1:
            pass
        else:
             dict1[key]=colors[len(dict1)%len(colors)]
    return dict1
def get_otu_color(taxonomy_file,phylum_colors,feature_table):
    otu_lineage = {}
    taxon = pd.read_csv(taxonomy_file, sep='\t')
    taxon = taxon.set_index(taxon['Feature ID'])
    otutable = biom.load_table(feature_table)
    otutable = otutable.to_dataframe()
    for ele in otutable.index:
        lineage = taxon.loc[ele]['Taxon']
        lineage = lineage.split(';')
        try:
            otu_lineage[ele] = lineage[1]
        except:
            otu_lineage[ele] = lineage[0]
    otu_color={}
    for key in otu_lineage:
        phylum = otu_lineage[key]
        otu_color[key] = phylum_colors[phylum] 
    return otu_color

def get_x_coordinates(tree, dist=1):
    """Get the x coordinates.
    Args:
        tree: a tree with sample_dict, which can be generate from
            the gini index compute module.
        corr:
            a method which compute the correlation between two arrays.
    Return:
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
        """
    ycoords ={}
    for clade in tree.find_clades(order='level'):
        #TODO corr function
        ycoords[clade] = clade.corr_value
        #corr(list(clade.sample_dict.values()))
    return ycoords
def new_get_lines(tree,xcoords,ycoords):
    traces = []
    for clade in tree.find_clades(order='level'):
        x0 = xcoords[clade]
        y0 = ycoords[clade]
        if hasattr(clade,'clades'):
            for ele in clade.clades:
                x_c = xcoords[ele]
                y_c = ycoords[ele]
                trace_p = go.Scatter(
                    x=(x0,x_c),
                    y=(y0,y_c),
                    marker=dict(color='rgb(0,0,0)',opacity=1),
                    showlegend=False,
                    mode = 'lines',
                    hoverinfo = 'none'
                )
                traces.append(trace_p)
    return traces

def get_lines(tree,clade,xcoords,ycoords,traces):
    path = tree.get_path(clade)
    x0 = xcoords[clade]
    y0 = ycoords[clade]
    if len(path)>1:
        parent= path[-2]
        x_p = xcoords[parent]
        y_p = ycoords[parent]
        trace_p = go.Scatter(# trace parent vertical line
            x = (x0,x_p),
            y = (y0,y_p),
            marker=dict(color='rgb(155,155,155)', opacity=0),
            showlegend=False,
            mode = 'lines',
            hoverinfo='none'
        )
        traces.append(trace_p)
    else:
        if len(path) == 1: # this calde is the subclade of root
            y_r = ycoords[tree.root]
            x_r = xcoords[tree.root]
            trace_p = go.Scatter(x=(x0,x_r),
                                y=(y0,y_r),
                                 mode='lines',
                                 marker =dict(color='rgb(125,125,125)',opacity=1),
                                 showlegend = False,
                                 hoverinfo='none'
                                )
            traces.append(trace_p)
            
    if clade.clades:
        '''
        x_l = xcoords[clade.clades[0]]
        x_r = xcoords[clade.clades[-1]]
        trace_c = go.Scatter(
        x = (x_l,x_r),
        y = (y0,y0) ,
        marker=dict(color='rgb(25,25,25)'),
        mode = 'lines',
        hoverinfo='none'
        )
        traces.append(trace_c)
        '''
        for ele in clade.clades:
            get_lines(tree,ele,xcoords,ycoords,traces=traces)
    else:
        pass
'''
        # dash lines
        if y0 >0:
            trace_t = go.Scatter(
            line = dict(color='rgb(25,25,25)',dash='dash'),
            mode = 'lines',
            hoverinfo='none',
            x =(x0,x0),
            y=(y0,0))
            traces.append(trace_t)
        '''
def obtain_series(metadata_file,col):
    df = pd.read_csv(metadata_file, sep='\t')
    df = df.set_index(df[df.columns[0]])
    series = df[col]
    return series

def compute_corr_value(tree, series,method='spearman'):
    #TODO judge the type of series value
    try:
        float(series[len(series)-1])
        is_num = True
    except:
        is_num = False
    corr_methods = {'spearman': scipy.stats.spearmanr,
                    'pearson': scipy.stats.pearsonr
    }
    stat_methods = {'t_test':stats_test.t_test,
                    'F_test':stats_test.F_test
    }
    corr = ''
    if is_num: #correlation
        if method in corr_methods:
            corr = corr_methods[method]
        else:
            corr = corr_methods['spearman']
    else:  # statistical method
        if method in stat_methods:
            corr = stat_methods[method]
        else:
            corr = stat_methods['t_test']
    for clade in tree.find_clades(order='level'):
        otu_list = list(clade.sample_dict.values())
        label_list = []
        for ele in clade.sample_dict:
            label_list.append(series[ele])
        if is_num:
            clade.corr_value = corr(otu_list,label_list)[0]
        else:
            clade.corr_value = corr(otu_list,label_list)
    return tree
def recurse_to_obtain_domain_otu(root_node):
    if not hasattr(root_node,'child_name'):
        tmp_sum = {}
        
        for ele in root_node.clades:
            seq, child_name = recurse_to_obtain_domain_otu(ele)
            tmp_sum[child_name] = seq #.seq
        root_node.child_name = max(tmp_sum, key=tmp_sum.get)
        root_node.seq = tmp_sum[root_node.child_name]
    return root_node.seq,root_node.child_name
def obtain_domain_otu(tree):
    """ Obtain the domain child for every node in the tree.
    (e.g, tree is (A,B)C;  A,B are OTUs, and A possess 8 seq , B 6 seq,
    we do C.child_name=A )
    Args:
        tree: a tree object which has been compute the sample_dict(by 
        gini_index_compute.py)
    """
    for clade in tree.get_terminals():
        clade.seq =sum(clade.sample_dict.values())
        clade.child_name = clade.name
    tmp = recurse_to_obtain_domain_otu(tree.root)
    return tree
    
def plot_tree(tree,otu_color, phylum_colors):
    text_dict ={}
    i = 0
    for clade in tree.find_clades(order='level'):
        if clade.name:
            text_dict[clade]='<br>name:'+clade.name+'<br> clade num:'+str(i)
        else:
            text_dict[clade]='<br>clade num:'+str(i)
        if clade.seq:
            text_dict[clade] += '<br>domian OTU seq: '+str(clade.seq)
        i+=1
    #tree = Phylo.read('test_tree.nwk','newick')
    y_coords = get_y_coordinates(tree)
    x_coords = get_x_coordinates(tree)
    X = []
    Y = []
    text = []
    colors = []
    for key in x_coords.keys():
        X.append(x_coords[key])
        Y.append(y_coords[key])
        text.append(text_dict[key])
        try:
            colors.append(otu_color[key.name])
        except:
            colors.append(otu_color[key.child_name])
    traces = new_get_lines(tree,x_coords,y_coords)
    #get_lines(tree,tree.root,x_coords,y_coords,traces=traces)
    trace = go.Scatter(
        x = X,
        y = Y,
        mode = 'markers',
        marker = dict(color=colors),
        showlegend = False,
        text= text,
        name=''
    )
    data = [trace]
    phlylum_legends = []
    for ele in phylum_colors:
        legend_trace = go.Scatter(
            x=[np.mean(X)],
            y=[np.mean(Y)],
            mode = 'markers',
            marker = dict(opacity=1, color= phylum_colors[ele]),
            name = ele
        )
        phlylum_legends.append(legend_trace)
    for ele in phlylum_legends:
        data.append(ele)
    layout = dict(showlegend=True,
            xaxis=dict(
                autorange=True,
                 showgrid=False,
                zeroline=False,
                showline=False,
                ticks='',
                showticklabels=False
                ),
            yaxis=dict(
                autorange=True,
                showgrid=False,
                zeroline=False,
                showline=False,
                ticks='',
                showticklabels=False
                ),
            hovermode = 'closest'
    )
    for ele in traces:
        data.append(ele)
    fig =go.Figure(data=data,layout=layout)
    tree_div = plot(fig,output_type='div')
    return tree_div
def run_this_script(tree, feature_table,metadata,col,taxonomy):
    #tree = gini_index_compute.get_tree(tree_file)
    tree = gini_index_compute.perform(tree, feature_table)
    series = obtain_series(metadata,col)
    tree = compute_corr_value(tree,series)
    tree = obtain_domain_otu(tree)
    phylum_colors = get_phylum_colors(taxonomy)
    otu_color = get_otu_color(taxonomy,phylum_colors,feature_table)
    div = plot_tree(tree,otu_color,phylum_colors)
    return div
if __name__ == "__main__":
    #tree = Phylo.read('test_tree.nwk','newick')
    tree = gini_index_compute.perform()
    f = open('tree_test.html','w')
    f.write(plot_tree(tree))