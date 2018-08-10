import pandas as pd
from Bio import Phylo
def binary_search(A,t):
    l=0
    r=len(A)-1
    A=sorted(A)
    while l<=r:
        mid = int((l+r)/2)
        #print(mid)
        if t>A[mid]:
            l=mid+1
        elif t<A[mid]:
            r=mid-1
        else:
            return mid
    return 'not_found'
def map_taxo_nwk(taxo_file,nwk_file):

    taxo =pd.read_csv(taxo_file,sep='\t')
    tree = Phylo.read(nwk_file,"newick")
    terminals = tree.get_terminals()
    taxo = taxo.sort_values('Feature ID')
    feature_id_col = taxo['Feature ID']
    for ele in terminals:
        found_index = binary_search(feature_id_col,ele.name)
        if found_index == 'not_found':
            pass
        else:
            temp =  taxo['Taxon'][found_index]
            temp = temp.replace('; ','#')
            ele.name = temp
    Phylo.write(tree,'mapped_tree.nwk','newick')
if __name__ == "__main__" :
    map_taxo_nwk('taxonomy.tsv','tree.nwk')




    