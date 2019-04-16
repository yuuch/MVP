from Bio import Phylo
tree = Phylo.read('tree.nwk','newick')



print(tree)

ts = tree.get_terminals()
i =0 
for ele in ts:
    ele.seq = i
    ele.child_name = ele.name
    i+=1

#for ele in ts

def recurse(root_node):
    if not hasattr(root_node,'child_name'):
        tmp_sum = {}
        
        for ele in root_node.clades:
            seq, child_name = recurse(ele)
            tmp_sum[child_name] = seq #.seq
            #tmp_sum[ele.child_name] = recurse(root_node)#ele.seq
        root_node.child_name = max(tmp_sum)
        root_node.seq = tmp_sum[root_node.child_name]
        print(tmp_sum)
        #print(root_node.name,root_node.seq)
    #print(root_node.name)
    return root_node.seq,root_node.child_name
recurse(tree.root)
print(tree)
