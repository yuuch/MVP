def read_metadata(filename,sep='\t'):
    """
    return a list of tuples.like:[(feature0,feature0_type),(feature1,feature1_type)]
    """
    if filename[-4:]=='.tsv':
        sep = '\t'
    elif filename[-4:]=='.csv':
        sep = ','
    f = open(filename)
    features = f.readline()
    features = features.split(sep)
    features[-1]=features[-1][:-1]
    feat_type = f.readline()
    feat_type = feat_type.split(sep)
    feat_type[-1]=feat_type[-1][:-1]
    result = list(zip(features,feat_type))[1:]
    return result
if __name__ == "__main__":
    print(read_metadata('demo_metadata.tsv'))
