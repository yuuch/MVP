    function send_heatmap(){
                var dict1 = {};
                dict1.feature0 = document.getElementById('mySelect0').value;
                dict1.feature1 = document.getElementById('mySelect1').value;
                dict1.feature2 = document.getElementById('mySelect2').value;
                dict1.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
                dict1.feature_table = 'MVP/upload_files/'+document.getElementById("feature_table").value;
                dict1.prevalence = document.getElementById('prevalence').value;
                dict1.abundance = document.getElementById('abundance').value;
                dict1.variance = document.getElementById('variance').value;
                //alert(typeof(JSON.stringify(dict1.metadata)));
                console.log(JSON.stringify(dict1))

            $.ajax(
                {
                    type: "POST",
                    url: "/graph/heat_map",
                    data: JSON.stringify(dict1),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        //alert("ajax");
                        //var result = $('<div />').append(data).find("#result").html();
                        //$('#result').style="width:1000px;height:800px";
                        $('#heatmap_result').html(data[0]);
                        //document.getElementById('result').innerHTML=data;
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
            )
            }

function read_metadata() {
              //features=['age','country','gender'];
              var features;
              var metadata_filename = document.getElementById('metadata').value;
              var demo_filename = 'MVP/upload_files/demo/demo_metadata.tsv';
              if (metadata_filename == ''){
                  metadata_filename =demo_filename;
              }else{
                  metadata_filename='MVP/upload_files/'+metadata_filename;
              }
              $.ajax({
                        type: "POST",
                        url: "/graph/return_string",
                        data: JSON.stringify(metadata_filename),
                        dataType: "json",
                        contentType: "application/json;charset=utf-8",
                        success: function(data){
                                      features=data;
                                      //alert(data[0]);
                                      var data_length = Object.keys(data).length;
                                      for(i=0;i<data_length;i++){
              var option = document.createElement("option");
                 temp = features[i];
              option.text = temp;
              x0.add(option);
              };
                                      },

                        failure: function(errMsg){alert(errMsg);}
                    }).responsetext;
              var x0 = document.getElementById("mySelect0");
              
          };

            function plot_tree(){
                var tree_paras={};
                tree_paras.tree_file = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
                tree_paras.node_num = document.getElementById('node_num').value;
                tree_paras.file_type = document.getElementById('tree_format').value;
                tree_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
                tree_paras.taxonomy_file = 'MVP/upload_files/'+document.getElementById('taxonomy_file').value;
                tree_paras.tree_type = document.getElementById('tree_type').value;
                tree_paras.feature0 = document.getElementById('mySelect0').value;
                tree_paras.feature1 = document.getElementById('mySelect1').value;
                tree_paras.feature2 = document.getElementById('mySelect2').value;
                tree_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
                tree_paras.feature_table = 'MVP/upload_files/'+document.getElementById("feature_table").value;
                console.log(JSON.stringify(tree_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_tree",
                    data: JSON.stringify(tree_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        //alert("ajax");
                        //var result = $('<div />').append(data).find("#result").html();
                        //$('#result').style="width:1000px;height:800px";
                        $('#tree_result').html(data[0]);
                        $('#anno_result').html(data[1]);
                        $('#heatmap_result').html(data[2]);
                        //document.getElementById('result').innerHTML=data;
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
                )
            }
        function plot_abun(){
            var abun_paras = {};
            abun_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            abun_paras.taxonomy_file = 'MVP/upload_files/'+document.getElementById('taxonomy_file').value;
            abun_paras.feature0 = document.getElementById('mySelect0').value;
            abun_paras.feature1 = document.getElementById('mySelect1').value;
            abun_paras.feature2 = document.getElementById('mySelect2').value;
            abun_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            abun_paras.feature_table = 'MVP/upload_files/'+document.getElementById("feature_table").value;
            abun_paras.abun_type = document.getElementById('abundance').value;
            abun_paras.log_flag = document.getElementById('log_flag').value;
            console.log(JSON.stringify(abun_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_abun",
                    data: JSON.stringify(abun_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        //alert("ajax");
                        //var result = $('<div />').append(data).find("#result").html();
                        //$('#result').style="width:1000px;height:800px";
                        $('#abun_result').html(data[0]);
                        $('#anno_result').html(data[1]);
                        $('#heatmap_result').html(data[2]);
                        //document.getElementById('result').innerHTML=data;
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
                )
            }
        function plot_stats_result(){
            var stats_paras = {};
            stats_paras.stats_method = document.getElementById('statistics methods').value;
            stats_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            stats_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            stats_paras.label_col = document.getElementById('mySelect0').value;
            console.log(JSON.stringify(stats_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_stats_test",// TODO
                    data: JSON.stringify(stats_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#stat_test_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
                )
        }

        function plot_osea_result(){
            var osea_paras = {};
            osea_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            osea_paras.taxonomy = 'MVP/upload_files/'+document.getElementById('taxonomy_file').value;
            osea_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            osea_paras.obj_col = document.getElementById('mySelect0').value;
            osea_paras.set_level = document.getElementById('set_level').value;
            console.log(JSON.stringify(osea_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_OSEA",// TODO
                    data: JSON.stringify(osea_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#osea_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
                )
        }
        function plot_dim_reduce_result(){
            var dim_paras = {};
            dim_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            dim_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            dim_paras.obj_col = document.getElementById('mySelect0').value;
            dim_paras.flag_3d = document.getElementById('flag_3d').value;
            dim_paras.n_component = document.getElementById('n_component').value;
            dim_paras.method = document.getElementById('method').value;
            console.log(JSON.stringify(dim_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_dim_reduce",// TODO
                    data: JSON.stringify(dim_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#dim_reduce_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems')
                    }
                }
                )
        }
               
        function plot_alpha_div(){
            var alpha_paras = {};
            alpha_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            alpha_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            alpha_paras.obj_col = document.getElementById('mySelect0').value;
            alpha_paras.tree = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
            alpha_paras.metric = document.getElementById('alpha_metric').value;
            console.log(JSON.stringify(alpha_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_alpha_diversity",// TODO
                    data: JSON.stringify(alpha_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#alpha_diversity_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems in plot alpha diversity')
                    }
                }
                )

        }
        function plot_beta_div(){
            var beta_paras = {};
            beta_paras.metric = document.getElementById('beta_metric').value;
            beta_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            beta_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            beta_paras.tree = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
            beta_paras.obj_col = document.getElementById('mySelect0').value;
            beta_paras.beta_dim_method = document.getElementById('beta_dim_method').value;
            beta_paras.n_components = document.getElementById('beta_n_components').value;
            console.log(JSON.stringify(beta_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_beta_diversity",// TODO
                    data: JSON.stringify(beta_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#beta_diversity_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems in plot beta diversity')
                    }
                }
                )

        }
        function plot_phylo_tree(){
            var beta_paras = {};
            beta_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            beta_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            beta_paras.obj_col = document.getElementById('mySelect0').value;
            beta_paras.node_num = document.getElementById('ecology_tree_node_num').value;
            beta_paras.tree = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
            beta_paras.taxonomy = 'MVP/upload_files/'+document.getElementById('taxonomy_file').value;
            console.log(JSON.stringify(beta_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_corr_tree",// TODO
                    data: JSON.stringify(beta_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#phylo_tree_result').html(data[0]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems in plot corr tree')
                    }
                }
                )



        }
        function plot_alpha_rarefaction(){
            var rare_paras = {};
            rare_paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            rare_paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            rare_paras.obj_col = document.getElementById('mySelect0').value;
            rare_paras.tree = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
            rare_paras.metric = document.getElementById('rarefaction_metric').value;
            rare_paras.step = document.getElementById('step').value;
            rare_paras.max_seq = document.getElementById('max_seq').value;
            rare_paras.rarefied_num = document.getElementById('rarefied_num').value;
            console.log(JSON.stringify(rare_paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_alpha_rarefaction",// TODO
                    data: JSON.stringify(rare_paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#rarefaction_box_result').html(data[0]);
                        $('#rarefaction_scatter_result').html(data[1]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems in alpha_rarefaction')
                    }
                }
                )

        }

function plot_ecology_scatter(){
            var paras = {};
            paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
            paras.feature_table = 'MVP/upload_files/'+document.getElementById('feature_table').value;
            paras.obj_col = document.getElementById('mySelect0').value;
            paras.tree = 'MVP/upload_files/'+document.getElementById('tree_file_name').value;
            paras.stats_method = document.getElementById('stats_test').value;
            paras.corr_method = document.getElementById('corr_coef').value;
            paras.ID_num = document.getElementById('ID_num').value;
            paras.taxonomy = 'MVP/upload_files/'+document.getElementById('taxonomy_file').value;
            console.log(JSON.stringify(paras))
                $.ajax({
                    type: "POST",
                    url: "/graph/plot_ecology_scatters",// TODO
                    data: JSON.stringify(paras),
                    dataType: "json",
                    contentType: "applications/json;charset=utf-8",
                    success: function(data){
                        $('#plot_subtree').html(data[0]);
                        $('#corr_GI_result').html(data[1]);
                        $('#corr_abu_result').html(data[2]);
                        $('#whole_abu_GI_result').html(data[3]);
                        $('#sub_tree_anno').html(data[4]);
                    },
                    error: function(xhr,status){
                        alert('sorry, there are problems in plot ecology scatter')
                    }
                }
                )
        }
function jump_html(){
    $.ajax({
        type: "POST",
        url: "/graph/jump_html"
    }
    )
}
function plot_PCA(){
        var paras = {};
        paras.metadata = 'MVP/upload_files/'+document.getElementById('metadata').value;
    $.ajax({
        type: "POST",
        url: "/graph/plot_PCA",
        data: JSON.stringify(paras),
        dataType: "json",
        contentType: "applications/json;charset=utf-8",
        success: function(data){
                     $('#PCA_result').html(data[0]);
        },
        error: function(xhr,status){
                    alert('sorry, there are problems in plot PCA')
                    }

    })
}