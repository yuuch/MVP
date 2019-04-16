def hex2rgb(hex_str):
    if len(hex_str) == 7:
        hex_str = hex_str[1:]
    
    r = int(hex_str[:2],16)
    g = int(hex_str[2:4],16)
    b = int(hex_str[4:6],16)
    return "'"+'rgb('+str(r)+','+str(g)+','+str(b)+')'+ "'"
if __name__ == "__main__":
    strs = '#42d4f4,#e6194b,#f58231,#3cb44b,#4363d8,#911eb4,#f032e6,#1f77b4,#ff7f0e,#2ca02c,#d62728,#9467bd,#8c564b,#e377c2,#7f7f7f,#bcbd22,#17becf' 
    strs = strs.split(',')
    result = ''
    for ele in strs:
        result += hex2rgb(ele)
        result +=','
    print(result)