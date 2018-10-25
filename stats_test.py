import numpy as np 
import scipy.stats as stats
import pandas as pd
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
    return method(arr1,arr2)

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