    #!/usr/bin/env python
    # coding: utf-8



from Bio import Phylo
import copy
from sklearn.utils import shuffle
import numpy as np
from sklearn import svm
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import cross_validate
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
import itertools
import pandas as pd
import biom
import time
import pickle
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
import corr_tree_new
import pandas as pd 
import copy 
import numpy as np 
import sklearn
import re

# def functions
def get_terminal_length(tree):
    terminal_lengths = {}
    for terminal in tree.get_terminals():
        path = tree.get_path(terminal)
        TBL = 0 # terminal branch length
        for ele in path:
            TBL += ele.branch_length
        terminal_lengths[terminal] = TBL
    sorted_terminal = sorted(terminal_lengths,key=terminal_lengths.get,reverse=True)
    sorted_terminal_lengths ={}
    for ele in sorted_terminal:
        sorted_terminal_lengths[ele] = terminal_lengths[ele]
   # print(sorted_terminal_lengths)
    return sorted_terminal_lengths

def cut_tree(tree,cp=0.08):
    """ 
    cut some branches of the tree.
    return the pruned branches.
    Args:
        tree: a Phylo.tree object
        cp: a percent used to cut the edge.
    """
    cp = 1-cp
    clades = []
    const = 1 # limit the terminal of a subtree
    terminal_lengths = get_terminal_length(tree)
    #print('terminal len:',len(terminal_lens))
    for terminal in terminal_lengths:
        current_terminals = tree.get_terminals()
        if not terminal in current_terminals:
            continue
        threshold = terminal_lengths[terminal]*cp
        temp = 0
        path = tree.get_path(terminal)
        for ele in path:
            temp += ele.branch_length
            if temp > threshold:
                #if not ele in current_terminals:  
                if len(ele.get_terminals()) > const:
                    backup = copy.copy(ele)
                    ele.clades = []
                    clades.append(backup)
                    try:
                        ele_parent = tree.get_path(ele)[-2]
                        ele_parent.clades.remove(ele)
                        if not ele_parent.clades:
                            ele_parent.sample_series *= 0
                    except:
                        print('need some other way to remove this clade')
                        ele.sample_series *= 0
                break
    return clades
def recurse_to_update(node):
    """used in the update sample series function.
    """
    if node.clades:
        node.sample_series = copy.deepcopy(recurse_to_update(node.clades[0]))
        for i,ele in enumerate(node.clades):
            if i == 0:
                continue
            else:
                node.sample_series += recurse_to_update(ele)
    else:
        return node.sample_series
def recursion_tree(node):
        """recursion to get the sample_series of a tree
        """
        if node.clades: # for non-leaf node
            tmp = 0
            flag = 0
            for clade in node.clades:
                if flag == 0:
                    tmp = copy.deepcopy(recursion_tree(clade).sample_series)
                else:
                    tmp += recursion_tree(clade).sample_series   
                flag = 1
            node.sample_series = copy.deepcopy(tmp)
        else: # leaf node which has been init above.
            try:
                a = node.sample_series
                #print(node.name +' is a leaf')
            except:
                print('please initialize the tree leaves by otu table.')
        return node    
def update_sample_series(tree):
    """ update sample series after cut the tree.
    Args:
        tree: tree with sample_siries on the terminals.
    """
    return recursion_tree(tree)
    #return tree
    
def get_mvp_tree(feature_table_path, tree_path, metadata_path):
    mvp_tree = corr_tree_new.MvpTree(feature_table_path, tree_path, metadata_path)
    return mvp_tree
def sep_mvp_tree(mvp_tree,cp):
    tree = copy.copy(mvp_tree.feature_tree)
    sub_trees = cut_tree(tree,cp)
    for ele in sub_trees:
        ele = update_sample_series(ele)
    print('len sub_trees:')
    print(len(sub_trees))
    return sub_trees

def get_node_score(node,coef=[1,1,1]):
    score = coef[0]*np.mean(node.sample_series)+        coef[1]*node.GI+coef[2]*np.std(node.sample_series)
    return score     

def search_sub_tree(sub_tree,coef,list1 = []):
    #print('search_sub_tree')
    p_score  = get_node_score(sub_tree,coef)# parent score
    child_score = 0
    for ele in sub_tree.clades:
        child_score += get_node_score(ele,coef)
    try:
        child_score = child_score/len(sub_tree.clades)
    except:
        pass
    if p_score <= child_score: 
        for ele in sub_tree.clades:
            search_sub_tree(ele,coef,list1)
    else:
        list1.append(sub_tree)
    #print('len list1')
    #print(len(list1))
    return list1

    
def generate_new_OTU_table(sub_trees,coef):
    """get a OTU table from the root node of the subtrees.
    Args:
        list1: obtained from the search_sub_tree function.
        coef: used for compute the score function,e.g.[1,1,1]
    """
    series = []
    for sub_tree in sub_trees:
        list1 = search_sub_tree(sub_tree,coef,list1=[])
        for ele in list1:
            series.append(ele.sample_series)
    df = pd.DataFrame(series)
    print('df.shape',df.shape)
    #print('df.index.name',df.index)
    return df


# useless functions?

def remove_non_otu_leaf(tree):
    """ for singular tree whose leaves contain some non otu node ,we need to delete the non otu node.
        or we can just change its sample_series to all zeros.
    """
    for t in tree.get_terminals():
        if t.name == None:
            t.sample_series *=0
def check_updated(subtree):
    flag = 0
    children = 0
    for ele in subtree.clades:
        if flag == 0:
            children = ele.sample_series
            flag = 1
        else:
            children+= ele.sample_series
    if subtree.sample_series.all() == children.all():
        try:
            for ele in subtree.clades:
                check_updated(ele)
        except:
            pass
    else:
        print('bad tree:')

def run_this_script(pickle_file, metadata_path, \
         obj_col, pos_label):

    # load mvp tree
    try:
        with open(pickle_file, 'rb') as f:  
            mvp_tree = pickle.load(f)
    except:
        mvp_tree = 0 #TODO
    backup_mvp_tree = copy.deepcopy(mvp_tree.feature_tree)
    # no  optmized ML
    b_auc, b_acc =no_optimized_ML(metadata_path,mvp_tree,obj_col,pos_label)
    # optimize and ML
    Max_key,best_result =cut_tree_and_ML(mvp_tree,metadata_path,  \
        obj_col, pos_label,backup_mvp_tree)
    a_auc = best_result['mean_auc']
    a_acc = best_result['mean_acc']
    mvp_tree.feature_tree = copy.deepcopy(backup_mvp_tree)
    df = obtain_optimized_df(Max_key,mvp_tree)
    return b_auc, b_acc, a_auc, a_acc,df

def no_optimized_ML(metadata_path,mvp_tree,obj_col,pos_label=None):
    SampleID = 'sample_name'
    metadata = pd.read_csv(metadata_path, sep= '\t')
    cols = metadata.columns
    p1 = '.*[Ss][Aa][Mm][Pp][Ll].*[IiNn][DdAa].*'
    pattern = re.compile(p1)
    for ele in cols:
        if len(pattern.findall(ele)) > 0:
            if len(pattern.findall(ele)[0])>5:
                SampleID = ele
    metadata = metadata.set_index(SampleID)
    if pos_label == None:
        pos_label = metadata[obj_col][0]
    # prepare data for ML
    feature_table_X = []
    feature_table_Y = []
    for ele in mvp_tree.feature_table.T.index:
        feature_table_X.append(mvp_tree.feature_table.T.loc[ele])
        if metadata[obj_col][ele] == pos_label:
            feature_table_Y.append(metadata[obj_col][ele])
        else:
            feature_table_Y.append('neg_label')
    feature_table_X = np.array(feature_table_X)
    feature_table_Y = np.array(feature_table_Y)
    skf1 = StratifiedKFold(n_splits=10,random_state = 1)
    skf1.get_n_splits(feature_table_X,feature_table_Y)

    old_clf = RandomForestClassifier(n_estimators=30, random_state=1)
    aucs = []
    fprs = []
    tprs = []
    uncut_accs = []
    for train_index, test_index in skf1.split(feature_table_X,feature_table_Y):
        X_train, X_test = feature_table_X[train_index], feature_table_X[test_index]
        y_train, y_test = feature_table_Y[train_index], feature_table_Y[test_index]
        old_clf.fit(X_train,y_train)
        y_pred_prob = old_clf.predict_proba(X_test)[:,0]

        y_test = [[0,1][ele == pos_label] for ele in y_test]
        y_pred = [[0,1][ele == pos_label] for ele in old_clf.predict(X_test)]

        auc = sklearn.metrics.roc_auc_score(y_test,y_pred_prob)
        fpr, tpr, _ = roc_curve(y_test,y_pred_prob)
        aucs.append(auc)
        fprs.append(fpr)
        tprs.append(tpr)
        acc = sklearn.metrics.accuracy_score(y_pred,y_test)
        uncut_accs.append(acc)
    mean_aucs = np.mean(aucs)
    mean_accs = np.mean(uncut_accs)
    return mean_aucs, mean_accs











def cut_tree_and_ML(mvp_tree,metadata_path,obj_col,\
        pos_label,backup_mvp_tree):
    skf = StratifiedKFold(n_splits=10,random_state=1)
    init_terminals_num = len(backup_mvp_tree.get_terminals())
    metadata = pd.read_csv(metadata_path, sep= '\t')
    SampleID = 'sample_name'
    metadata = pd.read_csv(metadata_path, sep= '\t')
    cols = metadata.columns
    p1 = '.*[Ss][Aa][Mm][Pp][Ll].*[IiNn][DdAa].*'
    pattern = re.compile(p1)
    for ele in cols:
        if len(pattern.findall(ele)) > 0:
            if len(pattern.findall(ele)[0])>5:
                SampleID = ele
    metadata = metadata.set_index(SampleID)
    tic = time.time()
    cps = [0.08]
    coefs = [
        [2,6,13]
        ]
    acc_threshold = 0.5            
    continue_flag = True
    cp_delta = 0.01
    deltas = [1,2,3]
    iter_num = 0
    # up down bounds used to constrain the cp(cut percent) value.
    up_bound = 0.15
    down_bound = 0.01
    # Max_key used to store the key.e.g str([0.08,1,20,1])
    Max_key = ''
    Last_Max_key = 'temp'
    results = {}
    #clfs = {}
    Max_df = 'placeholder'
    while continue_flag:
        iter_num += 1
        for cp in cps:
            mvp_tree.feature_tree =copy.deepcopy(backup_mvp_tree)
            if len(mvp_tree.feature_tree.get_terminals()) < init_terminals_num:
                print('tree has been edited')
                break
            sub_trees = sep_mvp_tree(mvp_tree,cp=cp)
            for coef in coefs:
                key = str([cp]+coef)
                if str([cp]+coef) in results:
                    continue
                else:
                    results[key] = {}
                df = generate_new_OTU_table(sub_trees,coef).T
                X = []
                Y = []
                if df.shape == (0,0):
                    continue
                for ele in df.index:
                    X.append(df.loc[ele])
                    if metadata[obj_col][ele] == pos_label:
                        Y.append(metadata[obj_col][ele])
                    else:
                        Y.append('neg_label')
                X = np.array(X)
                y = np.array(Y)
                clf = RandomForestClassifier(n_estimators=100, random_state=1)
                cut_aucs = []
                cut_fprs = []
                cut_tprs = []
                accs = []
                for train_index, test_index in skf.split(X,y):
                    X_train, X_test = X[train_index], X[test_index]
                    y_train, y_test = y[train_index], y[test_index]
                    clf.fit(X_train,y_train)
                    y_pred_prob = clf.predict_proba(X_test)[:,0]
                    y_test = [[0,1][ele==pos_label] for ele in y_test]
                    auc = sklearn.metrics.roc_auc_score(y_test,y_pred_prob)
                    y_pred = [[0,1][ele==pos_label] for ele in clf.predict(X_test)]
                    fpr, tpr, _ = roc_curve(y_test,y_pred_prob)
                    cut_aucs.append(auc)
                    cut_fprs.append(fpr)
                    cut_tprs.append(tpr)
                    assert len(y_test) ==len(y_pred)
                    acc = sklearn.metrics.accuracy_score(y_pred,y_test)
                    accs.append(acc)
                    results[key]['mean_auc'] = np.mean(cut_aucs)
                    results[key]['mean_acc'] = np.mean(accs)
                #print('mean auc:',np.mean(cut_aucs))
                #print('mean acc:',np.mean(accs))
        # generate new cps and coefs
        for key in results:
            if np.mean(results[key]['mean_acc']) > acc_threshold:
                acc_threshold = np.mean(results[key]['mean_acc'])
                Max_key = key
            else:
                pass
        if Max_key == Last_Max_key:
            print(Max_key)
            break
        Last_Max_key = copy.deepcopy(Max_key)
        Max_key_list = [float(ele) for ele in Max_key[1:-1].split(', ')]
        # update cps
        cps  = [max(Max_key_list[0]-cp_delta,down_bound), Max_key_list[0], min(Max_key_list[0]+cp_delta,up_bound)]
        # update coef
        Max_key_list.remove(Max_key_list[0])
        temps = [] # store paied coefs
        for i,ele  in enumerate(Max_key_list):
            temps.append([ele+deltas[i],ele-deltas[i]])
        coefs = list(itertools.product(temps[0],temps[1],temps[2]))
        coefs = [list(ele) for ele in coefs]
        #coef.append(temp)
        #print(results)
        #print('cps',cps)
        #print('coefs',coefs)
        #print('iter_num:',iter_num)
        new_tic = time.time()
        print('#seconds:',new_tic-tic)
        return Max_key, results[Max_key]

def obtain_optimized_df(Max_key,mvp_tree):
    """ Max key e.g str([0.08,1,2,5])
    """
    Max_key = [float(ele) for ele in Max_key[1:-1].split(', ')]
    cp = Max_key[0]
    coef = Max_key[1:]
    sub_trees = sep_mvp_tree(mvp_tree,cp=cp)
    df = generate_new_OTU_table(sub_trees,coef).T
    return df 

def get_index(y):
    pos_index = []
    neg_index = []
    for i,ele in enumerate(y):
        if ele == 'healthy':
            pos_index.append(i)
        else:
            neg_index.append(i)
    return pos_index + neg_index

if  __name__ == "__main__":

    pickle_file = 'pickles/dog_idb_tree.pickle'
    #feature_table_path = 'Dog_IDB/all.biom'
    tree_path = 'Dog_IDB/tree.nwk'
    metadata_path = 'Dog_IDB/833_20180418-110621.txt'
    obj_col = 'disease_stat'
    SampleID = 'sample_name'

    #mvp_tree = corr_tree_new.MvpTree(feature_table_path,tree_path,metadata_path)
    #with open(pickle_file,'wb') as f:
    #    pickle.dump(mvp_tree,f)
    results = run_this_script(pickle_file,metadata_path,\
        obj_col,pos_label='healthy')
    print(results)