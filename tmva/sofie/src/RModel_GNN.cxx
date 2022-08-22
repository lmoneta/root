#include <limits>
#include <algorithm>
#include <cctype>

#include "TMVA/RModel_GNN.hxx"


namespace TMVA{
namespace Experimental{
namespace SOFIE{


    RModel_GNN::RModel_GNN(RModel_GNN&& other){
      edges_block = std::move(other.edges_block);
      nodes_block = std::move(other.nodes_block);
      globals_block = std::move(other.globals_block);
      edges = std::move(other.edges);
      nodes = std::move(other.nodes);
      globals = std::move(other.globals);
      senders = std::move(other.senders);
      receivers = std::move(other.receivers);
   }

   RModel_GNN& RModel_GNN::operator=(RModel_GNN&& other){
      edges_block = std::move(other.edges_block);
      nodes_block = std::move(other.nodes_block);
      globals_block = std::move(other.globals_block);
      edges = std::move(other.edges);
      nodes = std::move(other.nodes);
      globals = std::move(other.globals);
      senders = std::move(other.senders);
      receivers = std::move(other.receivers);
   }

    RModel_GNN::RModel_GNN(const GNN_Input& graph_input){
        edges_block = std::make_unique<RFunction>(graph_input.edges_block);
        nodes_block = std::make_unique<RFunction>(graph_input.nodes_block);
        globals_block = std::make_unique<RFunction>(graph_input.globals_block);
        nodes = std::move(graph_input.nodes);
        globals = std::move(graph_input.globals);
        nodes = std::move(other.nodes);
        for(auto& it:graph_input.edges){
            senders.emplace_back(it.first);
            receivers.emplace_back(it.second);
        }
   }

    void RModel_GNN::AddFunction(std::unique_ptr<ROperator> func){
        switch(func->fTarget){
            case(FunctionTarget::NODES){
                nodes_block.emplace_back(func);
                break;
            }
            case(FunctionTarget::EDGES){
                edges_block.emplace_back(func);
                break;
            }
            case(FunctionTarget::GLOBALS){
                globals_block.emplace_back(func);
                break;
            }
        }
   }




}//SOFIE
}//Experimental
}//TMVA
