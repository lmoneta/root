#ifndef TMVA_SOFIE_RMODEL_GNN
#define TMVA_SOFIE_RMODEL_GNN

#include "TMVA/RModel.hxx"
#include "TMVA/RFunction.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

struct GNN_Input {
    RFunction edges_block;
    RFunction nodes_block;
    RFunction globals_block;
    std::vector<std::string> nodes;
    std::vector<std::pair<int,int>> edges;
    std::vector<std::string> globals;
};

class RModel_GNN: public RModel{

private:
   
    std::unique_ptr<RFunction> edges_block;
    std::unique_ptr<RFunction> nodes_block;
    std::unique_ptr<RFunction> globals_block;

    std::vector<std::pair<int,int>> edges; // contains node indices
    std::vector<std::string> nodes;
    std::vector<std::string> globals;
    std::vector<int> senders;              // contains node indices
    std::vector<int> receivers;            // contains node indices

public:

   //explicit move ctor/assn
   RModel_GNN(RModel_GNN&& other);

   RModel_GNN& operator=(RModel_GNN&& other);

   //disallow copy
   RModel_GNN(const RModel_GNN& other) = delete;
   RModel_GNN& operator=(const RModel_GNN& other) = delete;

   RModel_GNN(const GNN_Input& graph_input);
   RModel_GNN(){}
   RModel_GNN(std::string name, std::string parsedtime);

   void AddFunction(std::unique_ptr<RFunction> func);
   
   void InitializeGNN(int batch_size=1);
   void Generate(std::underlying_type_t<Options> options, int batchSize = 1);
   void Generate(Options options = Options::kDefault, int batchSize = 1) {
      Generate(static_cast<std::underlying_type_t<Options>>(options), batchSize);
   }
   
   ~RModel_GNN(){}}
   ClassDef(RModel_GNN,1);
};

}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RMODEL_GNN