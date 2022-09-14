#include <limits>
#include <algorithm>
#include <cctype>

#include "TMVA/RModel_GNN.hxx"


namespace TMVA{
namespace Experimental{
namespace SOFIE{


    RModel_GNN::RModel_GNN(RModel_GNN&& other){
      edges_update_block = std::move(other.edges_update_block);
      nodes_update_block = std::move(other.nodes_update_block);
      globals_update_block = std::move(other.globals_update_block);

      edge_node_agg_block = std::move(other.edge_node_agg_block);
      edge_global_agg_block = std::move(other.edge_global_agg_block);
      node_global_agg_block = std::move(other.node_global_agg_block);

      edges = std::move(other.edges);
      nodes = std::move(other.nodes);
      globals = std::move(other.globals);
      senders = std::move(other.senders);
      receivers = std::move(other.receivers);
    }

   RModel_GNN& RModel_GNN::operator=(RModel_GNN&& other){
      edges_update_block = std::move(other.edges_update_block);
      nodes_update_block = std::move(other.nodes_update_block);
      globals_update_block = std::move(other.globals_update_block);

      edge_node_agg_block = std::move(other.edge_node_agg_block);
      edge_global_agg_block = std::move(other.edge_global_agg_block);
      node_global_agg_block = std::move(other.node_global_agg_block);

      edges = std::move(other.edges);
      nodes = std::move(other.nodes);
      globals = std::move(other.globals);
      senders = std::move(other.senders);
      receivers = std::move(other.receivers);
      return *this;
    }

    RModel_GNN::RModel_GNN(const GNN_Init& graph_input_struct){
        edges_update_block.reset((graph_input_struct.edges_update_block).get());
        nodes_update_block.reset((graph_input_struct.nodes_update_block).get());
        globals_update_block.reset((graph_input_struct.globals_update_block).get());

        edge_node_agg_block.reset((graph_input_struct.edge_node_agg_block).get());
        edge_global_agg_block.reset((graph_input_struct.edge_global_agg_block).get());
        node_global_agg_block.reset((graph_input_struct.node_global_agg_block).get());

        nodes = std::move(graph_input_struct.nodes);
        globals = std::move(graph_input_struct.globals);
        nodes = std::move(graph_input_struct.nodes);
        for(auto& it:graph_input_struct.edges){
            senders.emplace_back(it.first);
            receivers.emplace_back(it.second);
        }
    }

    void RModel_GNN::GenerateGNN(int batchSize){
        Generate(Options::kGNN, batchSize);
        
        // computing inplace on input graph
        fGC += "GNN::GNN_Data infer(GNN::GNN_Data input_graph){\n";
        
        fGC+=edges_update_block->GenerateModel("Edge_Update",{{num_edge_features},{num_node_features},{num_node_features},{num_global_features}});
        fGC+=nodes_update_block->GenerateModel("Node_Update",{{num_edge_features+num_node_features+num_node_features},{num_node_features},{num_global_features}});
        fGC+=globals_update_block->GenerateModel("Global_Update",{{num_edge_features+num_node_features+num_node_features},{num_node_features},{num_global_features}});
        
        std::vector<std::vector<std::size_t>> AggregateInputShapes;
        for(int i=0; i<num_edges;++i){
            AggregateInputShapes[i] = {num_edge_features,num_node_features,num_node_features};
        }
        fGC+=edge_node_agg_block->GenerateModel("Edge_Node_Aggregate",AggregateInputShapes);
        fGC+=edge_global_agg_block->GenerateModel("Edge_Global_Aggregate",AggregateInputShapes);

        AggregateInputShapes.clear();
        for(int i=0; i<num_nodes;++i){
            AggregateInputShapes[i] = {num_node_features};
        }
        fGC+=node_global_agg_block->GenerateModel("Node_Global_Aggregate",AggregateInputShapes);

        // computing updated edge attributes
        for(int k=0; k<num_edges; ++k){
            fGC+="std::copy(input_graph.edge_data.begin()+"+std::to_string(num_edge_features)+"*"+std::to_string(k)+",input_graph.edge_data.begin()+"+std::to_string(num_edge_features)+"*"+std::to_string(k)+"+"+num_edge_features+",";
            fGC+=edges_update_block->Generate({"input_graph.edge_data.data()"+std::to_string(k),"input_graph.node_data.data()+input_graph.edges["+std::to_string(k)+"].first","input_graph.node_data.data()+input_graph.edges["+std::to_string(k)+"].second","input_graph.global_data.data()"});
        }

        for(int i=0; i<num_nodes; ++i){
            std::vector<GNN::GNN_Agg> agg_data_per_node;
            for(int k=0; k<num_edges; ++k){
                if(edges[k].first == i)
                    agg_data_per_node.push_back({nodes[i],nodes[edges[k].second]});
                else if(edges[k].second == i)
                    agg_data_per_node.push_back({nodes[i],nodes[edges[k].first]});
                else  
                    continue;
            }
            fGC+=edge_node_agg_block->Generate("Edge_Node_Agg"+i,agg_data_per_node,"Edge_Node_Agg"+std::to_string(i),batchSize);                      // aggregating edge attributes per node
            
            fGC+="std::copy(input_graph.node_data.begin()+input_graph.num_node_features*"std::to_string(i)+",input_graph.node_data.begin()+input_graph.num_node_features*"+std::to_string(i)+"+input_graph.num_node_features,";
            fGC+=nodes_update_block->Generate("Node_"+std::to_string(i)+"_Update",{"Edge_Node_Agg"+std::to_string(i),nodes[i],globals},"Node_"+std::to_string(i),batchSize);    // computing updated node attributes 
        }

        std::vector<GNN::GNN_Agg> agg_data;
        for(int k=0; k<edges.size(); ++k){
                agg_data.push_back({edges[k],nodes[i],nodes[edges[k].second]});
        }
        fGC+=edge_global_agg_block->Generate("Edge Global Aggregate",agg_data,"Edge_Global_Agg",batchSize);     // aggregating edge attributes globally
        fGC+=node_global_agg_block->Generate("Node Global Aggregate",agg_data,"Node_Global_Agg",batchSize);     // aggregating node attributes globally
        fGC+="input_graph.global_data=";
        fGC+=globals_update_block->Generate("Global Update",{"Edge_Global_Agg","Node_Global_Agg",globals},globals,batchSize); // computing updated global attributes

        fGC+="\nreturn input_graph;\n}";
        if (fUseSession) {
            fGC += "};\n";
         }
         fGC += ("} //TMVA_SOFIE_" + fName + "\n");
         fGC += "\n#endif  // TMVA_SOFIE_" + fName + "\n";
    }

    void RModel_GNN::AddFunction(std::unique_ptr<RFunction> func){
        if(func->GetFunctionType() == FunctionType::UPDATE){
            switch(func->GetFunctionTarget()){
                case(FunctionTarget::NODES){
                    nodes_update_block.reset(func.get());
                    break;
                }
                case(FunctionTarget::EDGES){
                    edges_update_block.reset(func.get());
                    break;
                }
                case(FunctionTarget::GLOBALS){
                    globals_update_block.reset(func.get());
                    break;
                }
            }
        } else{
            switch(func->GetFunctionRelation()){
                case(FunctionRelation::NODES_GLOBALS): {
                    node_global_agg_block.reset(func.get());
                    break;
                }
                case(FunctionRelation::EDGES_GLOBALS): {
                    edge_global_agg_block.reset(func.get());
                    break;
                }
                case(FunctionRelation::EDGES_NODES): {
                    edge_node_agg_block.reset(func.get());
                    break;
                }
            }
        }
    }




}//SOFIE
}//Experimental
}//TMVA
