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
import pandas as pd
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
import pickle
import perform_osea
import dimension_reduce
import diversity
import PCA_plot_ywch
from corr_tree_new import MvpTree
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
    ID_num =int(content['node_num'])
    #prevalence = content['prevalence']
    #abundance = content['abundance']
    #variance = content ['variance']
    try:
        f =  open('MVP/pickles/'+metadata.split('/')[-1] + '_heatmap.pickle','rb')
        heatmap_instance = pickle.load(f)
        print('read heatmap from pickle')
        f.close()
    except:
        heatmap_instance = heatmap.Heatmap(metadata, feature_table)
        heatmap_instance.map()
        with open('MVP/pickles/'+metadata.split('/')[-1] + '_heatmap.pickle','wb') as g:
            pickle.dump(heatmap_instance,g)
            print('write heatmap to pickle')
    #heatmap_instance.filter(prevalence_threshold=prevalence,abundance_num=abundance,variance_num=variance)
    heatmap_instance.map()
    heatmap_instance.sort_by_features(features[0],features[1],features[2])
    try:
        f = open('MVP/pickles/'+metadata.split('/')[-1]+'_mvp_tree.pickle','rb')
        mvp_tree = pickle.load(f)
        print('read mvp_tree from pickle')
        mvp_tree.get_subtree(ID_num)
        cols = [ele.name for ele in mvp_tree.subtree.get_terminals()]
        f.close()
    except:
        string_ = 'there are no pickles to read.please try plot_tree button'
        result = {0: string_}
        return jsonify(result)
    heatmap_instance.obtain_numerical_matrix(cols)
    show_label = content['show_label']
    if show_label == 'show': # show metadata besides the heatmap or not
        show_label = True
    else:
        show_label = False
    heatmap_div = heatmap_instance.plotly_div(show_label)
    result = {0:heatmap_div}
    return jsonify(result)

@bp.route('/plot_tree',methods=('GET','POST'))
def plot_tree():
    content = request.get_json(force=True)
    tree = content['tree_file']
    #tree_type = content['tree_type']
    #file_type = content['file_type']
    ID_num = int(content['node_num'])
    feature_table = content['feature_table']
    taxo_file = content['taxonomy_file']
    metadata  = content['metadata']
    #tree = circular_tree.read_tree(tree_file,file_type)
    try:
        f = open(metadata.split('/')[-1]+'_mvp_tree.pickle','rb')
        mvp_tree = pickle.load(f)
        print('read mvp_tree from pickle')
        f.close()
    except:
        mvp_tree = corr_tree_new.MvpTree(feature_table,tree,metadata,taxo_file,ID_num)
        file_paras = {'feature_table': feature_table,
            'metadata': metadata,
            'taxonomy': taxo_file,
            'tree': tree
            }
        with open('MVP/pickles/files.pickle','wb') as f:
            pickle.dump(file_paras,f)
        with open('MVP/pickles/'+metadata.split('/')[-1]+'_mvp_tree.pickle', 'wb') as g:
            pickle.dump(mvp_tree,g)
            print('wirte mvp_tree to pickle')
    mvp_tree.get_subtree(ID_num)
    cols = [ele.name for ele in mvp_tree.subtree.get_terminals()]
    ann = annotation.Annotation(cols,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    mvp_tree.get_colors(ann.colors,ann.mapped_phylum_colors)
    tree_div = mvp_tree.plot_tree()
    # plot_anno
    ann = annotation.Annotation(cols,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    
    #plot_heatmap
    features = [content['feature0'], content['feature1'], content['feature2']]
    try:
        f =  open('MVP/pickles/'+metadata.split('/')[-1] + '_heatmap.pickle','rb')
        heatmap_instance = pickle.load(f)
        print('read heatmap from pickle')
        f.close()
    except:
        heatmap_instance = heatmap.Heatmap(metadata, feature_table)
        heatmap_instance.map()
        with open('MVP/pickles/'+metadata.split('/')[-1] + '_heatmap.pickle','wb') as g:
            pickle.dump(heatmap_instance,g)
            print('write heatmap to pickle')
    heatmap_instance.sort_by_features(features[0], features[1], features[2])
    heatmap_instance.obtain_numerical_matrix(cols)

    show_label = content['show_label']
    if show_label == 'show': # show metadata besides the heatmap or not
        show_label = True
    else:
        show_label = False
    heatmap_div = heatmap_instance.plotly_div(show_label)
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
    obj_col = content['obj_col']
    stats_method = content['stats_method']
    corr_method = content['corr_method']
    ID_num = int(content['ID_num'])
    try:
        with open('MVP/pickles/files.pickle','rb') as f:
            files = pickle.load(f)
        feature_table = files['feature_table']
        tree = files['tree']
        taxo_file = files['taxonomy']
        metadata = files['metadata']
    except:
        print('no files.pickle exist please go to main page get main view first')
        pass

    try:
        f = open('MVP/pickles/'+metadata.split('/')[-1]+'_mvp_tree.pickle','rb')
        mvp_tree = pickle.load(f)
        print('read mvp_tree from pickle')
        f.close()
    except:
        mvp_tree = corr_tree_new.MvpTree(feature_table,tree,metadata,taxo_file,ID_num)
        file_paras = {'feature_table': feature_table,
            'metadata': metadata,
            'taxonomy': taxo_file,
            'tree': tree
            }
    mvp_tree.get_subtree(ID_num)
    cols = [ele.name for ele in mvp_tree.subtree.get_terminals()]
    ann = annotation.Annotation(cols,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    mvp_tree.get_colors(ann.colors,ann.mapped_phylum_colors)
    tree_div = mvp_tree.plot_tree()
    scatter_whole_tree = mvp_tree.plot_whole_tree()
    scatter_div1 = ''
    scatter_div2 = ''
    if stats_method != 'None':
        mvp_tree.stats_test(obj_col,stats_method,ID_num)
        scatter_div1 = mvp_tree.plot_scatter('pvalue', 'GI',ID_num)
        scatter_div2 = mvp_tree.plot_scatter('pvalue', 'abundance',ID_num)
    if corr_method != 'None':
        mvp_tree.get_corr_coefficient(obj_col, corr_method,ID_num)
        scatter_div1 = mvp_tree.plot_scatter('corr_coef', 'GI',ID_num)
        scatter_div2 = mvp_tree.plot_scatter('corr_coef', 'abundance',ID_num)
    result = {0:tree_div,1:scatter_div1,2:scatter_div2,3:scatter_whole_tree,4:ann_div}
    return jsonify(result)

@bp.route('jump_html',methods=('GET','POST'))
def jump_html():
    return render_template('test_js.html')

@bp.route('/plot_PCA',methods=('GET','POST'))
def plot_PCA():
    content = request.get_json(force=True)
    print(content)
    metadata = content['metadata']
    tree_path = content['tree_path']
    feature_table = content['feature_table']
    ID_num =int(content['ID_num'])
    metadata_df = pd.read_csv(metadata,sep='\t')
    metadata_df = metadata_df.drop(0)
    #print(dir(PCA_plot_ywch))
    div = PCA_plot_ywch.run_this_script(metadata_df, feature_table, tree_path,\
        metadata, ID_num)
    result = {0:div}
    return jsonify(result)
@bp.route('/update_file_names', methods=['GET','POST'])
def update_file_names():
    try:
        content = request.get_json(force=True)
        print(content)
    except:
        pass
    files_dict = {}
    try:
        # metadata
        with open('MVP/pickles/metadata_filename.pickle','rb')as f:
            metadata_file = pickle.load(f)
            files_dict['metadata'] = metadata_file['metadata_filename']
        # tree
        with open('MVP/pickles/tree_filename.pickle','rb')as f:
            tree_file = pickle.load(f)
            files_dict['tree_file_name'] = tree_file['tree_filename']
        # taxonomy
        with open('MVP/pickles/taxonomy_filename.pickle','rb')as f:
            metadata_file = pickle.load(f)
            files_dict['taxonomy'] = metadata_file['taxonomy_filename']
            # TODO add taxonomy
        # feature_table
        with open('MVP/pickles/feature_table_filename.pickle', 'rb')as f:
            metadata_file = pickle.load(f)
            files_dict['feature_table'] = metadata_file['feature_table_filename']
    except:
        files_dict = {
            'feature_table':'feature-table.biom',
            'metadata': 'demo_metadata.tsv',
            'taxonomy': 'taxonomy.tsv',
            'tree_file_name':'tree.nwk'
        }
    return jsonify(files_dict)
@bp.route('/reload_metadata', methods=['GET', 'POST'])
def reload_metadata():
    try:
        with open('MVP/pickles/files.pickle','rb') as f:
            files = pickle.load(f)
            metadata = files['metadata']
    except:
        metadata = 'MVP/upload_files/demo_metadata.tsv'
    meta_data_list = read_metadata.read_metadata(metadata)
    d1 ={}
    for i in range(len(meta_data_list)):
        d1[i]=meta_data_list[i]
    jsd1= jsonify(d1)
    #print(jsd1)
    return jsd1
