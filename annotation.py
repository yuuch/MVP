import biom
import pandas as pd 
from Bio import Phylo 
import plotly
import plotly.graph_objs as go
import copy
class Annotation(object):
    def __init__(self, cols, feature_table_file, taxo_file):
        #self.tree = tree#Phylo.read(file = tree_file,format = tree_file_format)
        self.cols = cols
        self.feature_table = biom.load_table(feature_table_file).to_dataframe().to_dense()
        self.taxonomy = pd.read_csv(taxo_file,sep='\t')
        self.taxonomy.index = self.taxonomy['Feature ID']
        self.map_taxo(self.cols)
        self.generate_colors()

    def map_taxo(self,cols):
        #terminals = tree.get_terminals()
        #feature_names = [node.name for node in terminals]
        feature_names = cols
        feature_name_count_dict = {}
        for feature_name in feature_names:
            # feature_count: the number of seq in one feature(otu)
            #print(self.feature_table.loc[feature_name])
            feature_count = 0
            if len(self.feature_table.loc[feature_name])==0:
                continue
            else:
                feature_count = sum(list(self.feature_table.loc[feature_name]))
            unassigned = ['Unassigned',' p__unassigned',' c__unassigned',' o__unassigned',' f__unassigned'
                         ,' g__unassigned',' s__unassigned']
            if feature_count == 0:
                continue
            taxo_of_feature = self.taxonomy.loc[feature_name]['Taxon'].split(';')
            while len(taxo_of_feature)<len(unassigned):
                taxo_of_feature.append(unassigned[len(taxo_of_feature)])
            levels = {}
            for i in range(7):
                levels['level'+str(i)]={}
            for i in range(len(taxo_of_feature)):
                if taxo_of_feature[i] in levels['level'+str(i)]:
                    levels['level'+str(i)][taxo_of_feature[i]] += feature_count
                else:
                    levels['level'+str(i)][taxo_of_feature[i]] = feature_count
            
            #for i in range(len(levels)):
                #if levels['level'+str(i)] == {}:
                    #levels['level'+str(i)] = unassigned[i]
            #while len(levels) < len(unassigned):
                #levels['level'+str(len(levels))] = unassigned[len(levels)]    
            feature_name_count_dict[feature_name]=levels
        self.barplot_dict = feature_name_count_dict
    '''
    def add_node_num(self,tree):
        nodes = tree.get_terminals()+tree.get_non_terminals()
        i = 0
        for node in nodes:
            node.node_num = i 
            i += 1
    '''
    def generate_colors(self):
        #colors =[]
        colors = ['rgb(255,0,0)','rgb(255,247,0)','rgb(255,0,247)','rgb(162,255,0)',
        'rgb(255,111,0)','rgb(111,0,255)','rgb(2,104,23)','rgb(104,2,40)','rgb(2,23,104)',
        'rgb(192,77,11)','rgb(12,192,162)','rgb(126,3,145)','rgb(17,102,7)','rgb(57,12,179)',
        'rgb(195,216,8)','rgb(216,8,154)','rgb(6,97,94)','rgb(97,6,36)','rgb(125,218,5)',
        'rgb(255,0,80)']
        self.colors = colors
    def plot_annotation(self):
        mapped_level={'level0':'Kingdom','level1':'Phylum','level2':'Class','level3':'Order',
                 'level4':'Family','level5':'Genus','level6':'Species','level7':'OTU'}
        data = []
        level_appeard_name = {}
        color_index = 0
        for i in range(7):
            level_appeard_name['level'+str(i)]={}
        for otu in self.barplot_dict:
            otu_data = []
            #reverse_level=[]
            #for level in self.barplot_dict[otu]:
                #reverse_level.append(level)
            #reverse_level.reverse()
            for level in self.barplot_dict[otu]:
                for ele in self.barplot_dict[otu][level]:
                    temp_index = 0
                    not_appeared = True
                    if ele in level_appeard_name[level]:
                        temp_index = level_appeard_name[level][ele]
                        not_appeared = False
                    else:
                        #color_index +=1
                        #color_index %= len(self.colors)
                        temp_index = (len(level_appeard_name[level])+int(level[5])*3)%len(self.colors)
                        level_appeard_name[level][ele]=temp_index
                    temp_color = self.colors[temp_index]
                    if ele =='Unassigned' or ele[3:]=='unassigned':
                        temp_color = 'rgb(0,0,0)'
                    filt_legend=['level1','level2']
                    filt_result=False
                    #print(ele[0])
                    if level in filt_legend:
                        filt_result=True
                    trace=go.Bar(y=[mapped_level[level]],
                                     x = [self.barplot_dict[otu][level][ele]],
                                     orientation='h',
                                     name = ele,
                                     showlegend = not_appeared and filt_result,
                                     legendgroup = mapped_level[level],
                                     marker = dict(color=temp_color)
                        )
                    otu_data.append(trace)
            if len(data)>0:
                last_otu_data = copy.copy(data[-1])
                to_be_del=[]
                for i in range(7):
                    #print(min(len(otu_data),len(last_otu_data)))
                    #to_be_del = []
                    if i < min(len(otu_data),len(data[-1])):
                        if otu_data[i].name == data[-1][i].name:
                            to_be_del.append(i)
                            otu_data[i].x=tuple((data[-1][i].x[0]+otu_data[i].x[0],))
                            data[-1][i].x=tuple((0,))
            else:
                pass
            data.append(otu_data)
        #for i
        self.mapped_phylum_colors = level_appeard_name['level1']
        result_data = []
        for i in range(7):
            trace = go.Bar(y=[mapped_level['level'+str(6-i)]],
                           x=[0],
                           orientation='h',
                           name=mapped_level['level'+str(6-i)],
                           showlegend=False
                          )
            result_data.append(trace)
        #print(data[0])
        for i in range(7):
            for ele in data:
                if i<len(ele):
                    if ele[i].x[0]>0:
                        result_data.append(ele[i])
        #print(result_data)
        layout = go.Layout(
            #autosize =False,
            #height = 500,
            #width = 1000,
            barmode = 'stack',
            xaxis=dict(
                autorange=True,
                 showgrid=False,
                zeroline=False,
                showline=False,
                ticks='',
                showticklabels=False
                ),
            hovermode='closest'
        )
        fig = go.Figure(data=result_data,layout=layout)
        div_str = plotly.offline.plot(fig,filename='stack_bar.html',output_type='div')
        return div_str