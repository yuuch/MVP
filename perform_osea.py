import time
import numpy as np
import os
import OSEA
import pandas as pd
import biom
import heatmap
import stats_test
import plotly
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
def run_osea(taxonomy_file, feature_table,metadata_file,obj_col,set_level,per_num=1000):
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
    rank_list = obtain_rank_list(part1,part2)
    osea_real = OSEA.OSEA(rank_list,Taxon_file=taxonomy_file,set_level=set_level)
    ### permutation many times to generate the null distribution
    enrichment_scores = [] # 
    tmp_df = heatmap_instance.df[heatmap_instance.df_primary_col]
    for i in range(per_num):
        rank_list = OSEA.permutation_to_obtain_ranklist(tmp_df)
        es = osea_real.get_ES(rank_list)
        enrichment_scores.append(es)
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
def plot_final_result(final_result):
    """ plot the result of the osea.
    
    """
    data = [plotly.graph_objs.Bar(
        x = list(final_result.keys()),
        y = [ele['ES'] for ele in final_result.values()],
        text = ['pvalue: '+str(ele['pvalue']) for ele in final_result.values()]
    )]
    layout = plotly.graph_objs.Layout(
        title = 'OSEA result'
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