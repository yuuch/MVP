from Bio import Phylo
import biom
def get_tree(tree_file):
    tree = Phylo.read(tree_file,"newick")
    return tree

def series2dict(series):
    dict1 = {}
    for i in range(len(series)):
        dict1[series.index[i]]=series[i]
    return dict1
def get_otus(OTU_table):
    """ read the otu table.
        Args:
            OTUS_table: a biom format OTU table.
        Return:
            return a dict include otus and samples.like:{
                 OTU0:{Sample0:12,Sample1:22,Sample3:2},
                 OTU1:{Sample0:2,Sample1:22,Sample3:22},
                 OTU2:{Sample0:133,Sample1:122,Sample3:52},
            }
    """
    table = biom.load_table(OTU_table)
    df = table.to_dataframe().transpose().to_dense()
    otus = {}
    for otu in df.columns:
        tmp = df[otu]
        tmp = series2dict(tmp)
        otus[otu] = tmp
    return otus


def initialize_tree_leaves(tree,leaves_dict):
    terminals = tree.get_terminals()
    for leaf in terminals:
        #print(leaf.name)
        leaf.sample_dict = leaves_dict[leaf.name] # sample_number is a dict
        leaf.sequence_counter = 0 
        temp=0
        total = 0
        for ele in leaf.sample_dict:
            leaf.sequence_counter+=leaf.sample_dict[ele]
            temp += leaf.sample_dict[ele]**2
            total+= leaf.sample_dict[ele]
            leaf.gini_index = 1-(temp*1.0/total**2)




def check_sample_dict_exist_in_clades(clade):
    label = 1
    for element in clade.clades:
        if hasattr(element,'sample_dict'):
            pass
        else:
            label*=0
    return label


def merge_dict(big_dict,small_dict):
    for ele in small_dict:
        if ele in big_dict:
            big_dict[ele]+=small_dict[ele]
        else:
            big_dict[ele]=small_dict[ele]
    return big_dict



def generate_sample_dict(ywch_clade):
    sample_dict = {}
    for ele in ywch_clade.clades:
        sample_dict=merge_dict(sample_dict,ele.sample_dict)
    return sample_dict
    




def get_gini_index(sample_dict):
    """ compute gini index for single node."""
    total =0
    square_total = 0
    for ele in sample_dict:
        total += sample_dict[ele]
        square_total+=(sample_dict[ele]**2)
    gini_index = 1-(square_total*1.0/total**2)
    return gini_index
        

def compute_recursion(root_node):
    if not check_sample_dict_exist_in_clades(root_node):
        for ele in root_node.clades:
            if not hasattr(ele,'sample_dict'):
                compute_recursion(ele)
    root_node.sample_dict = generate_sample_dict(root_node)
    root_node.gini_index = get_gini_index(root_node.sample_dict)
    


def get_tree_with_gini_index(tree,leaves_dict):
    """Compute the gini index of every node of the tree.
    Args:
        tree: a phylotree read by biopython.phylo.read.
        leaves_dict: a dict include the information of every leaf.For single 
        leaf,it is {OTU_0:12,OTU_1:123,OTU_4:121}
    Return:
        retur a tree with gini index in its nodes.
    """
    initialize_tree_leaves(tree,leaves_dict)
    root_node = tree.clade
    compute_recursion(root_node)
    return tree
#if __name__ == "__main__":
def perform(tree_file, feature_table):
    tree = get_tree(tree_file)
    otus = get_otus(feature_table)
    tree = get_tree_with_gini_index(tree,otus)
    return tree
    #print(tree)
#tree =get_tree_with_gini_index(tree,leaves_dict)
#print(tree)