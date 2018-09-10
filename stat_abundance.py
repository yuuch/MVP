import plotly
import plotly.graph_objs as go
import biom
import math
def plot_stat_abun(feature_table='feature-table.biom', abundance_type='absolute',log_flag='yes'):
    result_dict = {}
    df = biom.load_table(feature_table).to_dataframe().to_dense().transpose()
    for col in df.columns:
        result_dict[col]=sum(df[col])
    # sorted by abundance
    keys_list = sorted(result_dict,key=result_dict.get,reverse=True)
    new_result = {}
    for key in keys_list:
        new_result[key] = result_dict[key]
    result_dict = new_result
    # choose relative abundance or absolute abundance
    total = sum(result_dict.values())
    if abundance_type == 'relative':
        for ele in result_dict:
            result_dict[ele] /=total
    # log or not
    if log_flag == 'yes':
        for ele in result_dict:
            result_dict[ele] = math.log2(result_dict[ele]+2)
    trace = go.Bar(
        x = list(result_dict.keys()),
        y = list(result_dict.values())
    )
    layout = go.Layout(title=abundance_type +' abundace bar plot')
    fig = go.Figure(data=[trace],layout=layout)
    div = plotly.offline.plot(fig,output_type='div')
    return div,result_dict
if __name__ == '__main__':
    div = plot_stat_abun()
    f = open('stat_abun.html','w')
    f.write(div[0])
    
    


