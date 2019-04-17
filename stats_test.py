import numpy as np 
import scipy.stats as stats
import pandas as pd
import plotly
import math
def jacknife(arr):
    n = len(arr)
    if n == 1:
        return -1
    tmp = []
    for i in range(n):
        tmp.append((sum(arr)-arr[i])/(n-1))
    return tmp
        
def F_test(arr1, arr2, alpha=0.05):
    """Do F_test for two arrays. 
    We play a trick to make sure the big var divide the small var so that
    we can do the one tail F-test.
    """
    m = len(arr1)
    n = len(arr2)
    F = np.var(arr1)/np.var(arr2)
    if F < 1:
        F = 1/F  # make sure that the big var divide the lsmall var
        m, n = n, m
    pvalue = 1-stats.f.cdf(F, m-1, n-1)
    return pvalue

def t_test(arr1, arr2):
    return stats.ttest_ind(arr1, arr2)[1]
def mannwhitneyu_test(arr1,arr2):
    return stats.mannwhitneyu(arr1,arr2)[1]

def correlation(arr1, arr2, corr_type='spearman'):
    """get the correlation betweent two  arrays which possess the same length

    Return:
        correlation value: emmm..
        pvalue: the null hypothesis is arr1 and arr2 are not uncorrelated.
        For example, it will return (0.9,0.0) which means thre corr is 0.9 and
        arr1 and arr2 are significantly correlated.
        
    """
    methods = {'spearman':stats.spearmanr,
               'pearson':stats.pearsonr
               }
    method = methods[corr_type]
    return method(arr1,arr2)[0]
def fisher_exact_test(arr11, arr12, arr21, arr22):
    tmp11 = sum(arr11)
    tmp12 = sum(arr12)
    tmp21 = sum(arr21)
    tmp22 = sum(arr22)
    oddsratio,pvalue = stats.fisher_exact([[tmp11, tmp12], [tmp21, tmp22]])
    return oddsratio, pvalue

def choose_two_class(dataframe, column):
    tmp = ''
    #condition=lambda x:True if x==expect_label else False
    part1_index = []
    part2_index = []
    i=0
    for ele in dataframe[column]:
        if tmp == '':
            tmp = ele
        if ele == tmp:
            part1_index.append(i)
        else:
            part2_index.append(i)
        i+=1
    df_part1 = dataframe.iloc[part1_index]
    df_part2 = dataframe.iloc[part2_index]
    return df_part1,df_part2
def perform_test(df1, df2, method_name):
    """ perform test...for multi otus
        Args:
            df1, df2: dataframe1, 2 which possess the same columns
            method_name: a test name(string) like t_test of F_test or  fisher_exact_test
    """
    methods = {'F_test':F_test,
               't_test':t_test,
               'fisher_exact_test':fisher_exact_test}
    if method_name not in methods:
        print('Please input the right method name:t_test,F_test,fisher_exact_test')
        return -1
    method = methods[method_name]
    result_dict = {}
    if method_name == 'fisher_exact_name':
        pass
    else: # for t_test and F_test
        for col in df1.columns:
            value = method(df1[col], df2[col])
            result_dict[col] = -math.log(value)
    return result_dict
def plot_result_dict(result_dict,cols=None,taxo_file=None,colors=None,color_index=None):
    taxonomy = pd.read_csv(taxo_file, sep='\t')
    print('read files successfully')
    assert len(colors)
    data = []
    taxonomy = pd.read_csv(taxo_file, sep='\t')
    taxo_col = 'Taxon' #TODO re match
    feature_id = 'Feature ID'
    taxonomy = taxonomy.set_index(feature_id)
    if not cols:
        cols = result_dict.keys()
    for ele in cols:
        lineage =taxonomy.loc[ele][taxo_col].split(';') 
        try:
            phylum = lineage[1]
        except:
            phylum = ' p__unassigned'
        try:
            color = colors[color_index[phylum]]  
        except:
            color = 'rgb(0,0,0)'
        trace = plotly.graph_objs.Scatter(
            x = [ele],
            y = [result_dict[ele]],
            mode = 'markers',
            marker = dict(color=color),
            showlegend=False
        )
        data.append(trace)
    line_001 = plotly.graph_objs.Scatter(
        x = list(result_dict.keys()),
        y = [-math.log(0.01) for ele in result_dict.keys()],
        mode = 'lines',
        line=dict(color='rgb(0,0,125)',dash='dot'),
        name = '- log(0.01)'
    )
    data.append(line_001)
    line_005 = plotly.graph_objs.Scatter(
        x = list(result_dict.keys()),
        y = [-math.log(0.05) for ele in result_dict.keys()],
        mode = 'lines',
        line=dict(color='rgb(0,125,0)',dash='dash'),
        name = '- log(0.05)'
    )
    data.append(line_005)
    
    layout =plotly.graph_objs.Layout(
        autosize = True,
        yaxis=dict(title= '- log pvalue')

    )
    fig = plotly.graph_objs.Figure(data=data, layout=layout)
    div_str = plotly.offline.plot(fig, output_type='div')
    return div_str   
def run_this_script(dataframe,label_column,method_name):
    df1, df2 = choose_two_class(dataframe,label_column)
    result_dict = perform_test(df1, df2, method_name)
    return result_dict


