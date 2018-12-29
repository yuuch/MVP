import functools
import os
import plotly
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for,jsonify,json
)
from werkzeug.security import check_password_hash, generate_password_hash
import sys
import re
#print(sys.path)
sys.path.append("MVP/")
import alpha_rarefaction
import read_metadata
import heatmap
import circular_tree
import corr_tree
import corr_tree_new
import annotation
import rectangle_tree
import stat_abundance  
import stats_test
import OSEA    
import perform_osea
import dimension_reduce
import diversity
#import corr_tree
#from MVP.db import get_db
#url_for('static', filename='base.css')
#url_for('static', filename='base.js')
bp = Blueprint('graph', __name__, url_prefix='/graph')
@bp.route('/graph',methods=('GET','POST'))
def graph():
    return render_template('graph/test.html')
#def generate_new_bar():
  #  path = os.getcwd()+"/graph/bar.html"
    #new_bar = open(path)
@bp.route('/return_string',methods=('GET','POST'))
def return_string():
    d = {'test0':10989,'test1':222}
    content = request.get_json()
    #cwd = os.getcwd()
    f = os.path.join('',content)
    print(f)
    meta_data_list = read_metadata.read_metadata(f)
    #print(meta_data_list)
    #print(content)
    #d = {'features':meta_data_list}
    #print(d)
    #print(l)
    d1 ={}
    for i in range(len(meta_data_list)):
        d1[i]=meta_data_list[i]
    jsd1= jsonify(d1)
    #print(jsd1)
    return jsd1
@bp.route('/heat_map',methods=('GET','POST'))
def heat_map():
    content = request.get_json(force=True) # content is filename string
    metadata = content['metadata']
    feature_table = content['feature_table']
    features = [content['feature0'],content['feature1'],content['feature2']]
    prevalence = content['prevalence']
    abundance = content['abundance']
    variance = content ['variance']
    heatmap_instance = heatmap.Heatmap(metadata,feature_table)
    #heatmap_instance.filter(prevalence_threshold=prevalence,abundance_num=abundance,variance_num=variance)
    heatmap_instance.map()
    heatmap_instance.sort_by_features(features[0],features[1],features[2])
    cols = list(heatmap_instance.df.columns)
    heatmap_instance.obtain_numerical_matrix(cols)
    result = {0:heatmap_instance.plotly_div()}
    return jsonify(result)
@bp.route('/plot_tree',methods=('GET','POST'))
def plot_tree():
    content = request.get_json(force=True)
    tree_file = content['tree_file']
    tree_type = content['tree_type']
    file_type = content['file_type']
    node_num = int(content['node_num'])
    feature_table = content['feature_table_file']
    taxo_file = content['taxonomy_file']
    tree = circular_tree.read_tree(tree_file,file_type)
    # plot_tree
    sub_tree = circular_tree.obtain_subtree(tree,node_num)
    cols = [ele.name for ele in sub_tree.get_terminals()]
    tree_div = ''
    if tree_type == 'circular_tree':
        tree_div = circular_tree.plot_tree(sub_tree)
    else:
        tree_div = rectangle_tree.plot_tree(sub_tree)
    # plot_anno
    ann = annotation.Annotation(cols,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    #plot_heatmap
    metadata = content['metadata']
    feature_table = content['feature_table']
    features = [content['feature0'],content['feature1'],content['feature2']]
    heatmap_instance = heatmap.Heatmap(metadata,feature_table)
    heatmap_instance.map()
    heatmap_instance.sort_by_features(features[0],features[1],features[2])
    heatmap_instance.obtain_numerical_matrix(cols)
    heatmap_div = heatmap_instance.plotly_div()
    # total         
    result = {0:tree_div,1:ann_div,2:heatmap_div}
    return jsonify(result)
@bp.route('/plot_abun',methods=('GET','POST'))
def plot_abun():
    content = request.get_json(force=True)
    feature_table = content['feature_table']
    log_flag =content['log_flag']
    abun_type = content['abun_type']
    # abun
    abun_div_and_dict  = stat_abundance.plot_stat_abun(feature_table,abun_type,log_flag)
    abun_div = abun_div_and_dict[0]
    cols = [ele for ele in abun_div_and_dict[1]]
    # heatmap
    metadata = content['metadata']
    feature_table = content['feature_table']
    features = [content['feature0'],content['feature1'],content['feature2']]
    heatmap_instance = heatmap.Heatmap(metadata,feature_table)
    heatmap_instance.map()
    heatmap_instance.sort_by_features(features[0],features[1],features[2])
    heatmap_instance.obtain_numerical_matrix(cols)
    heatmap_div = heatmap_instance.plotly_div()
    # annotation
    taxo_file = content['taxonomy_file']
    ann = annotation.Annotation(cols,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    result = {0:abun_div,1:ann_div,2:heatmap_div}
    return jsonify(result)
@bp.route('/plot_stats_test',methods=('GET','POST'))
def plot_stats_test():
    content = request.get_json(force=True)
    label_col = content['label_col']
    metadata = content['metadata']
    feature_table = content['feature_table']
    test_method = content['stats_method']
    #features = [content['feature0'],content['feature1'],content['feature2']]
    heatmap_instance = heatmap.Heatmap(metadata,feature_table)
    heatmap_instance.map()
    part1, part2 = stats_test.choose_two_class(heatmap_instance.df,label_col)
    part1  = part1[heatmap_instance.df_primary_col]
    part2  = part2[heatmap_instance.df_primary_col]
    test_result = stats_test.perform_test(part1,part2,test_method)
    div_str = stats_test.plot_result_dict(test_result)
    result={0:div_str}
    return jsonify(result)

@bp.route('/plot_OSEA',methods=('GET','POST'))
def plot_OSEA():
    content = request.get_json(force=True)
    metadata =content['metadata']
    taxonomy = content['taxonomy']
    feature_table = content['feature_table']
    set_level = content['set_level']
    obj_col = content['obj_col']
    osea_result = perform_osea.run_osea(taxonomy, feature_table, 
        metadata,obj_col,set_level)
    result = {0:perform_osea.plot_final_result(osea_result)}
    return jsonify(result)

@bp.route('/plot_dim_reduce',methods=('GET' ,'POST'))
def plot_dim_reduce():
    content = request.get_json(force=True)
    metadata =content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    # new buttons
    n_component = int(content['n_component'])
    method = content['method']
    flag_3d = content['flag_3d']

    #print(type(n_component))
    #print(n_component)

    heatmap_instance = heatmap.Heatmap(metadata,feature_table)
    heatmap_instance.map()
    labels = heatmap_instance.df[obj_col]
    matrix = heatmap_instance.df[heatmap_instance.df_primary_col]
    reduced = dimension_reduce.reduce_dimension(matrix,n_component,method)
    for ele in reduced:
        print(len(ele))
        break
    div = dimension_reduce.dimension_reduce_visualize(reduced,labels,flag_3d)
    result = {0:div}
    return jsonify(result)
@bp.route('plot_corr_tree', methods=('GET', 'POST'))
def plot_corr_tree():
    content = request.get_json(force=True)
    metadata =content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    tree_file = content['tree']
    node_num = int(content['node_num'])
    taxonomy = content['taxonomy']

    tree = circular_tree.read_tree(tree_file)#  ,file_type)
    tree = circular_tree.obtain_subtree(tree,node_num)
    div = corr_tree.run_this_script(tree,feature_table,metadata,obj_col,taxonomy)
    result = {0:div}
    return jsonify(result)

@bp.route('plot_alpha_diversity',methods=('GET', 'POST'))
def plot_alpha_diversity():
    content = request.get_json(force=True)
    metadata = content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    metric = content['metric']
    tree = content['tree']
    alpha_table = diversity.alpha_diversity_pre(feature_table,metric,tree)
    tmp_result = diversity.alpha_diversity(
        alpha_table=alpha_table,
        metadata=metadata,
        label_col=obj_col)
    div = diversity.alpha_box_plot(tmp_result)
    result = {0:div}
    return jsonify(result)

@bp.route('plot_beta_diversity',methods=('GET', 'POST'))
def plot_beta_diversity():
    content = request.get_json(force=True)
    metadata = content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    tree = content['tree']
    metric = content['metric']
    dim_method = content['beta_dim_method']
    n_components = int(content['n_components'])
    distance_matrix =diversity.beta_diversity_pre(feature_table,tree,metric)

    beta_dict, axis_names = diversity.beta_diversity(
        col=obj_col,
        metadata_file=metadata,
        distance_matrix=distance_matrix,
        dim_method=dim_method,
        n_components=n_components
    )
    div = diversity.plot_beta_scatter(beta_dict,axis_names)
    result = {0:div}
    return jsonify(result)

@bp.route('plot_alpha_rarefaction',methods=('GET', 'POST'))
def plot_alpha_rarefaction():
    content = request.get_json(force=True)
    metadata = content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    metric = content['metric']
    tree = content['tree']
    rarefied_num = int(content['rarefied_num'])
    max_seq = int(content['max_seq'])
    step = int(content['step'])
    box, scatter = alpha_rarefaction.plot_alpha_rarefaction(
        feature_table,metadata,max_seq,step,metric,obj_col,rarefied_num,tree)
    result = {0: box, 1: scatter}
    return jsonify(result)

@bp.route('plot_ecology_scatters',methods=('GET', 'POST'))
def plot_ecology_scatters():
    content = request.get_json(force=True)
    metadata = content['metadata']
    feature_table = content['feature_table']
    obj_col = content['obj_col']
    tree = content['tree']
    stats_method = content['stats_method']
    corr_method = content['corr_method']
    ID_num = int(content['ID_num'])
    mvp_tree = corr_tree_new.MvpTree(feature_table,tree,metadata,ID_num)
    scatter_div1 = ''
    scatter_div2 = ''
    if stats_method != 'None':
        mvp_tree.stats_test(obj_col,stats_method)
        scatter_div1 = mvp_tree.plot_scatter('pvalue', 'GI')
        scatter_div2 = mvp_tree.plot_scatter('pvalue', 'abundance')
    if corr_method != 'None':
        mvp_tree.get_corr_coefficient(obj_col, corr_method)
        scatter_div1 = mvp_tree.plot_scatter('corr_coef', 'GI')
        scatter_div2 = mvp_tree.plot_scatter('corr_coef', 'abundance')
    tree_div = mvp_tree.plot_tree()
    scatter_whole_tree = mvp_tree.plot_whole_tree()
    result = {0:tree_div,1:scatter_div1,2:scatter_div2,3:scatter_whole_tree}
    return jsonify(result)
