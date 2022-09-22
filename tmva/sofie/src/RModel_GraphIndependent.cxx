#include <limits>
#include <algorithm>
#include <cctype>

#include "TMVA/RModel_GraphIndependent.hxx"


namespace TMVA{
namespace Experimental{
namespace SOFIE{


    RModel_GraphIndependent::RModel_GraphIndependent(RModel_GraphIndependent&& other){
      edges_update_block = std::move(other.edges_update_block);
      nodes_update_block = std::move(other.nodes_update_block);
      globals_update_block = std::move(other.globals_update_block);

      num_nodes = std::move(other.num_nodes);
      num_edges = std::move(other.num_edges);

      fName = std::move(other.fName);
      fFileName = std::move(other.fFileName);
      fParseTime = std::move(other.fParseTime);
    }

   RModel_GraphIndependent& RModel_GraphIndependent::operator=(RModel_GraphIndependent&& other){
      edges_update_block = std::move(other.edges_update_block);
      nodes_update_block = std::move(other.nodes_update_block);
      globals_update_block = std::move(other.globals_update_block);

      num_nodes = std::move(other.num_nodes);
      num_edges = std::move(other.num_edges);

      fName = std::move(other.fName);
      fFileName = std::move(other.fFileName);
      fParseTime = std::move(other.fParseTime);

      return *this;
    }

    RModel_GraphIndependent::RModel_GraphIndependent(const GraphIndependent_Init& graph_input_struct){
        edges_update_block.reset((graph_input_struct.edges_update_block).get());
        nodes_update_block.reset((graph_input_struct.nodes_update_block).get());
        globals_update_block.reset((graph_input_struct.globals_update_block).get());

        num_nodes = graph_input_struct.num_nodes;
        num_edges = graph_input_struct.edges.size();
        num_node_features = graph_input_struct.num_node_features;
        num_edge_features = graph_input_struct.num_edge_features;
        num_global_features = graph_input_struct.num_global_features;

        fFileName = graph_input_struct.filename;
        fName = fFileName.substr(0, fFileName.rfind("."));

        std::time_t ttime = std::time(0);
        std::tm* gmt_time = std::gmtime(&ttime);
        fParseTime  = std::asctime(gmt_time);
    }

    void RModel_GraphIndependent::Generate(int batchSize){
        std::string hgname;
        GenerateHeaderInfo(hgname);

        //Generating Infer function definition for Edge update function
        long next_pos;
        fGC+="\n\nnamespace Edge_Update{\n";
        std::vector<std::vector<std::size_t>> Update_Input = {{num_edge_features,1}};
        edges_update_block->Initialize();
        edges_update_block->AddInputTensors(Update_Input);
        fGC+=edges_update_block->GenerateModel(fName);
        next_pos = edges_update_block->GetFunctionBlock()->WriteInitializedTensorsToFile(fName);
        fGC+="}\n";

        fGC+="\n\nnamespace Node_Update{\n";
        // Generating Infer function definition for Node Update function
        Update_Input = {{num_node_features,1}};
        nodes_update_block->Initialize();
        nodes_update_block->AddInputTensors(Update_Input);
        fGC+=nodes_update_block->GenerateModel(fName,next_pos);
        next_pos = nodes_update_block->GetFunctionBlock()->WriteInitializedTensorsToFile(fName);
        fGC+="}\n";

        fGC+="\n\nnamespace Global_Update{\n";
        // Generating Infer function definition for Global Update function
        Update_Input = {{num_global_features,1}};
        globals_update_block->Initialize();
        globals_update_block->AddInputTensors(Update_Input);
        fGC+=globals_update_block->GenerateModel(fName,next_pos);
        next_pos = globals_update_block->GetFunctionBlock()->WriteInitializedTensorsToFile(fName);
        fGC+="}\n";
        
        // computing inplace on input graph
        fGC += "GNN::GNN_Data infer(GNN::GNN_Data input_graph){\n";
        
        // computing updated edge attributes
        for(int k=0; k<num_edges; ++k){
            fGC+="std::vector<float> Edge_"+std::to_string(k)+"_Update = ";
            fGC+=edges_update_block->Generate({"input_graph.edge_data.data()+"+std::to_string(k)});
            fGC+="\nstd::copy(Edge_"+std::to_string(k)+"_Update.begin(),Edge_"+std::to_string(k)+"_Update.end(),input_graph.edge_data.begin()+"+std::to_string(k)+");";
        }
        fGC+="\n";

        // computing updated node attributes
        for(int k=0; k<num_nodes; ++k){
            fGC+="std::vector<float> Node_"+std::to_string(k)+"_Update = ";
            fGC+=nodes_update_block->Generate({"input_graph.node_data.data()+"+std::to_string(k)});
            fGC+="\nstd::copy(Node_"+std::to_string(k)+"_Update.begin(),Node_"+std::to_string(k)+"_Update.end(),input_graph.node_data.begin()+"+std::to_string(k)+");";
        }
        fGC+="\n";

        // computing updated global attributes
        fGC+="input_graph.global_data=";
        fGC+=globals_update_block->Generate({"input_graph.global_data"}); 
        fGC+="\n";
        
        fGC+="\nreturn input_graph;\n}";
        fGC += ("} //TMVA_SOFIE_" + fName + "\n");
        fGC += "\n#endif  // TMVA_SOFIE_" + hgname + "\n";

    }

}//SOFIE
}//Experimental
}//TMVA