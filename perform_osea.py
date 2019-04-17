import time
import numpy as np
import os
import OSEA
import pandas as pd
import math
import biom
import heatmap
import stats_test
import plotly
from multiprocessing import Pool
from OSEA import permutation_to_obtain_ranklist
def obtain_rank_list(part1, part2, test_method_name='t_test'):
    methods = {'t_test': stats_test.t_test,
               'F_test': stats_test.F_test}
    test_method = methods[test_method_name]
    rank_list_unsort = {}
    for col in part1.columns:
        tmp_pvalue = test_method(part1[col],part2[col])
        rank_list_unsort[col]=tmp_pvalue
    tmp_list = sorted(rank_list_unsort, key=rank_list_unsort.get,reverse=True)
    rank_list = {}
    for ele in tmp_list:
        rank_list[ele]=rank_list_unsort[ele]
    return rank_list
def run_osea(taxonomy_file, feature_table,metadata_file,obj_col,set_level, \
        test_method_name='t_test',per_num=1000):
    """Generate an osea instance ,and compute it's Enrichment scores and pvalue.
    Args:
        obj_col:object column,which can be seen as a 0-1 label,dividing the
            samples into two groups.
        set_level: it it a string "Phylum", "Genus", and so on .
        per_num :permutation number
    Return:
        final_result: it is a dict object like {set:{ES: , pvalue: }}
    """
    heatmap_instance = heatmap.Heatmap(metadata_file,feature_table)
    heatmap_instance.map()
    part1, part2 = stats_test.choose_two_class(heatmap_instance.df,obj_col)
    part1 = part1[heatmap_instance.df_primary_col]
    part2 = part2[heatmap_instance.df_primary_col]
    rank_list = obtain_rank_list(part1,part2,test_method_name)
    osea_real = OSEA.OSEA(rank_list,Taxon_file=taxonomy_file,set_level=set_level)
    ### permutation many times to generate the null distribution
    enrichment_scores = [] # 
    tmp_df = heatmap_instance.df[heatmap_instance.df_primary_col]
    # TODO we can use pool here to accerate the program.(done)
    # multiprocess by using Pool()
    with Pool() as p:
        rank_lists = p.starmap(permutation_to_obtain_ranklist, \
            [[tmp_df, i, test_method_name] for i in range(per_num)])
    for rank_list in rank_lists:
        es = osea_real.get_ES(rank_list)
        enrichment_scores.append(es)
    """ single process code
    for i in range(per_num):
        rank_list = OSEA.permutation_to_obtain_ranklist(tmp_df,test_method_name)
        es = osea_real.get_ES(rank_list)
        enrichment_scores.append(es)
    """
    ### get the null distribution for every set
    set_enrichment_score = {}
    for key in enrichment_scores[0]:
        set_enrichment_score[key] = []
        for i in range(per_num):
            set_enrichment_score[key].append(enrichment_scores[i][key])
    ### compute the p-value
    final_result = {}
    for ele in set_enrichment_score:
        distribution = OSEA.generate_distribution(set_enrichment_score[ele])
        pvalue = OSEA.p_value(osea_real.es[ele],distribution)
        final_result[ele]={'ES':osea_real.es[ele],'pvalue':pvalue}
    return final_result
def get_phylumn_info(lineages,taxo_info):
    """get the phylumn info if you know the below level taxo info such as genus specis and so.on.
    Args:
        lineages: all OTUs' taxonomy infomation 
    """
    for lineage in lineages:
        lineage = lineage.split(';')
        if taxo_info in lineage and len(lineage) > 1:
            return lineage[1]
    return ' p__unassigned'

def plot_final_result(final_result, taxo_file,colors= None, color_index=None):
    """ plot the result of the osea.
    """
    taxo_all = pd.read_csv(taxo_file,sep='\t')
    columns = taxo_all.columns
    # TODO use re to match 'taxon' 'taxonomy' and so on
    lineages = taxo_all['Taxon']
    data = []
    try:
        assert not colors
        trace_bar = plotly.graph_objs.Bar(
            x = list(final_result.keys()),
            y = [-math.log(ele['pvalue'])for ele in final_result.values()],
            text = ['ES: '+str(ele['ES']) for ele in final_result.values()]
        )
        data.append(trace_bar)
    except:
        for ele in final_result.keys():
            ele_phylumn = get_phylumn_info(lineages,ele)
            value = final_result[ele]
            temp_color = 'rgb(0,0.0)'
            if ele_phylumn in color_index:
                temp_color = colors[color_index[ele_phylumn]]
            temp_trace = plotly.graph_objs.Bar(
                x=[ele],
                y=[-math.log(value['pvalue'])],
                text = ['ES:'+str(value['ES'])],
                marker=dict(color=temp_color),
                showlegend=False
                )
            data.append(temp_trace)
        print('there are colors')
    line_001 = plotly.graph_objs.Scatter(
        x = list(final_result.keys()),
        y = [-math.log(0.01) for ele in final_result.keys()],
        mode = 'lines',
        line=dict(color='rgb(0,0,125)',dash='dot'),
        name = '- log(0.01)'
    )
    data.append(line_001)
    line_005 = plotly.graph_objs.Scatter(
        x = list(final_result.keys()),
        y = [-math.log(0.05) for ele in final_result.keys()],
        mode = 'lines',
        line=dict(color='rgb(0,125,0)',dash='dash'),
        name = '- log(0.05)'
    )
    data.append(line_005)
    """
    trace_bar = plotly.graph_objs.Bar(
        x = list(final_result.keys()),
        y = [-math.log(ele['pvalue'])for ele in final_result.values()],
        text = ['ES: '+str(ele['ES']) for ele in final_result.values()]
    )
    """

    #########trace_scatter = plotly.graph_objs.Scatter(mode='line',)
    # plot line for -math.log(0.05 or 0.0.1)
    layout = plotly.graph_objs.Layout(
        title = 'OSEA result',
        yaxis = dict(title='- log pvalue')

    )
    fig = plotly.graph_objs.Figure(data=data, layout=layout)
    div_str = plotly.offline.plot(fig, output_type='div')
    return div_str

if __name__ == "__main__":
    taxonomy_file = os.path.join('','upload_files/taxonomy.tsv')
    feature_table = os.path.join('','upload_files/feature-table.biom')
    metadata_file = os.path.join('', 'upload_files/demo/demo_metadata.tsv')
    final_result = run_osea(taxonomy_file, feature_table, metadata_file,
                            obj_col='Subject', set_level='Genus')
    print(final_result)