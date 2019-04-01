import pandas as pd
import sklearn
from sklearn import svm
import numpy as np
from sklearn.utils import shuffle
def bfs(tree,threshold):
    latest_terminal = [tree.root]
    while True:
        temp = []
        for terminal in latest_terminal:
            new_terminal = []
            try:
                for clade in terminal.clades:
                    if score(clade)>threshold:
                        temp.append(clade)
                        new_termianl.append(clade)
            except:
                print('no clades')
            terminal = new_terminal
        if not temp:
            break
        else:
            print(temp)
            latest_terminal = temp
    return tree
bfs(tree)
def get_score(node,coefs):
    score = coefs[0]*np.mean(node.sample_series)+coefs[1]*node.GI+coefs[2](node.sample_series)
    return score
def tree_optimize(mvp_tree,coefs=None):
    """ prune a mvp tree.
    """
    if not coefs:
        coefs = [1,1,1]
        # TODO
    

    


def tree2OTU_table(mvp_tree):
    """ get an OTU table from a mvp_tree
    """
    series = []
    for terminal in mvp_tree.feature_tree.get_terminals():
        try:
            series.append(terminal.sample_series)
        except:
            print('there is no sample series in tree2OTU ')
    df = pd.dataframe(series)
    return df

def OTU_table_ML(OTU_table,metadata,obj_col):
    """ using ML technich to predict obj_col.
    args:
        OTU_table: a dataframe
        metadata: a dataframe
        obj_col: a string in the metadata columns.
    """
    for ele in OTU_table.index:
        #print(ele)
        X.append(df.loc[ele])
        Y.append(metadata[obj_col][ele])
    precisions = []
    for train_time in range(100):    
        X,Y = shuffle(X,Y)
        sample_num = len(X)
        sep_num = int(0.8*sample_num)
        train_set = [X[:sep_num],Y[:sep_num]]
        test_set = [X[sep_num:],Y[sep_num:]]
        clf = svm.SVC(gamma='scale')
        clf.fit(train_set[0], train_set[1])  
        predict_result = clf.predict(test_set[0])
        count = 0
        for i in range(len(predict_result)):
            if predict_result[i] == test_set[1][i]:
                count += 1
            else:
                pass
        precisions.append(1.0*count/len(predict_result))
    print(np.mean(precisions))

    

if __name__ == "__main__":
    coefs = [1,1,1]
    mvp_tree ='mvp_Tree'#TODO
    mvp_tree_pruned = tree_optimize(mvp_tree)
    OTU_table = tree2OTU_table(mvp_tree)
    precision = OTU_table_ML(OTU_table)


