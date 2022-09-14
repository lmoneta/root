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

      num_nodes = std::move(other.num_nodes);
      num_edges = std::move(other.num_edges);
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

      num_nodes = std::move(other.num_nodes);
      num_edges = std::move(other.num_edges);
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

        num_nodes = std::move(graph_input_struct.num_nodes);
        num_edges = std::move(graph_input_struct.edges.size());
        for(auto& it:graph_input_struct.edges){
            senders.emplace_back(it.first);
            receivers.emplace_back(it.second);
        }
    }

    void RModel_GNN::GenerateGNN(int batchSize){
        Generate(Options::kGNN, batchSize);
        
        // computing inplace on input graph
        fGC += "GNN::GNN_Data infer(GNN::GNN_Data input_graph){\n";
        
        // Generating Infer function definition for Edge Update function
        std::vector<std::vector<std::size_t>> Update_Input = {{num_edge_features},{num_node_features},{num_node_features},{num_global_features}};
        fGC+=edges_update_block->GenerateModel("Edge_Update", Update_Input);

        // Generating Infer function definition for Node Update function
        Update_Input = {{num_edge_features+num_node_features+num_node_features},{num_node_features},{num_global_features}};
        fGC+=nodes_update_block->GenerateModel("Node_Update", Update_Input);

        // Generating Infer function definition for Global Update function
        Update_Input = {{num_edge_features+num_node_features+num_node_features},{num_node_features},{num_global_features}};
        fGC+=globals_update_block->GenerateModel("Global_Update",Update_Input);
        
        std::vector<std::vector<std::vector<std::size_t>>> AggregateInputShapes;
        for(int i=0; i<num_edges;++i){
            AggregateInputShapes[i] = {{num_edge_features},{num_node_features},{num_node_features}};
        }
        fGC+=edge_node_agg_block->GenerateModel("Edge_Node_Aggregate",AggregateInputShapes);
        fGC+=edge_global_agg_block->GenerateModel("Edge_Global_Aggregate",AggregateInputShapes);

        AggregateInputShapes.clear();
        for(int i=0; i<num_nodes;++i){
            AggregateInputShapes[i] = {{num_node_features}};
        }
        fGC+=node_global_agg_block->GenerateModel("Node_Global_Aggregate",AggregateInputShapes);

        // computing updated edge attributes
        for(int k=0; k<num_edges; ++k){
            fGC+="std::vector<float> Edge_"+std::to_string(k)+"Update = ";
            fGC+=edges_update_block->Generate({"input_graph.edge_data.data()"+std::to_string(k),"input_graph.node_data.data()+input_graph.edges["+std::to_string(k)+"].first","input_graph.node_data.data()+input_graph.edges["+std::to_string(k)+"].second","input_graph.global_data.data()"});
            fGC+="std::copy(Edge_"+std::to_string(k)+"Update.begin(),Edge_"+std::to_string(k)+"Update.end(),input_graph.edge_data.begin()"+std::to_string(k)+");";
        }

        std::vector<std::string> Node_Edge_Aggregate_String;
        for(int i=0; i<num_nodes; ++i){
            for(int k=0; k<num_edges; ++k){
                if(receivers[k] == i){
                    Node_Edge_Aggregate_String.emplace_back("input_graph.edge_data.data()+"+std::to_string(k));
                    Node_Edge_Aggregate_String.emplace_back("input_graph.node_data.data()+"+std::to_string(i));
                    Node_Edge_Aggregate_String.emplace_back("input_graph.node_data.data()+"+std::to_string(senders[k]));
                }
            }
            fGC+="std::vector<float> Node_"+std::to_string(i)+"_Edge_Aggregate = ";
            fGC+=edge_node_agg_block->Generate(Node_Edge_Aggregate_String);                      // aggregating edge attributes per node

            fGC+="std::vector<float> Node"+std::to_string(i)+"_Update = ";
            fGC+=nodes_update_block->Generate({"Node_"+std::to_string(i)+"_Edge_Aggregate.data()","input_graph.node_data.data()"+std::to_string(i),"input_graph.global_data.data()"});    // computing updated node attributes 
            fGC+="std::copy(Node_"+std::to_string(i)+"_Edge_Update.begin(), Node_"+std::to_string(i)+"_Edge_Update.end(), input_graph.node_data.begin()+"+std::to_string(i);
            Node_Edge_Aggregate_String.clear();
        }

        std::vector<std::string> Node_Global_Aggregate_String;
        for(int k=0; k<num_nodes; ++k){
            Node_Global_Aggregate_String.emplace_back("input_graph.node_data.data()"+std::to_string(k));
        }

        std::vector<std::string> Edge_Global_Aggregate_String;
        for(int k=0; k<num_edges; ++k){
            Edge_Global_Aggregate_String.emplace_back("input_graph.edge_data.data()"+std::to_string(k));
            Edge_Global_Aggregate_String.emplace_back("input_graph.node_data.data()"+std::to_string(receivers[k]));
            Edge_Global_Aggregate_String.emplace_back("input_graph.node_data.data()"+std::to_string(senders[k]));
        }

        fGC+="std::vector<float> Edge_Global_Aggregate = ";
        fGC+=edge_global_agg_block->Generate(Edge_Global_Aggregate_String);     // aggregating edge attributes globally

        fGC+="std::vector<float> Node_Global_Aggregate = ";
        fGC+=node_global_agg_block->Generate(Node_Global_Aggregate_String);     // aggregating node attributes globally
        fGC+="input_graph.global_data=";
        fGC+=globals_update_block->Generate({"Edge_Global_Aggregate","Node_Global_Aggregate", "input_graph.global_data"}); // computing updated global attributes

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
                case(FunctionTarget::NODES): {
                    nodes_update_block.reset(func.get());
                    break;
                }
                case(FunctionTarget::EDGES): {
                    edges_update_block.reset(func.get());
                    break;
                }
                case(FunctionTarget::GLOBALS): {
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
