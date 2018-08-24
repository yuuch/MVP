import functools
import os
import plotly
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for,jsonify,json
)
from werkzeug.security import check_password_hash, generate_password_hash
import sys
#print(sys.path)
sys.path.append("MVP/")
import read_metadata
import heatmap
import circular_tree
import annotation
#from MVP.db import get_db

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
    heatmap_instance.filter(prevalence_threshold=prevalence,abundance_num=abundance,variance_num=variance)
    heatmap_instance.map()
    heatmap_instance.sort_by_features(features[0],features[1],features[2])
    heatmap_instance.obtain_numerical_matrix()
    result = heatmap_instance.plotly_div()
    return jsonify(result)
@bp.route('/plot_tree',methods=('GET','POST'))
def plot_tree():
    content = request.get_json(force=True)
    tree_file = content['tree_file']
    file_type = content['file_type']
    node_num = int(content['node_num'])
    feature_table = content['feature_table_file']
    taxo_file = content['taxonomy_file']
    tree = circular_tree.read_tree(tree_file,file_type)
    sub_tree = circular_tree.obtain_subtree(tree,node_num)
    tree_div = circular_tree.plot_tree(sub_tree)
    ann = annotation.Annotation(sub_tree,feature_table,taxo_file)
    ann_div = ann.plot_annotation()
    result = {0:tree_div,1:ann_div}
    return jsonify(result)

