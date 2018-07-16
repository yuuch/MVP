a = [1,2,3,5,8]
from time import time
f = open('matrix_test','w')
for ele in a:

    f.write(str(ele)+',')
    f.write(str(time())+"\n")
