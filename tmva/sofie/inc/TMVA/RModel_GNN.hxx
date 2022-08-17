#ifndef TMVA_SOFIE_RMODEL_GNN
#define TMVA_SOFIE_RMODEL_GNN

#include "TMVA/RModel.hxx"
#include "TMVA/RFunction.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel_GNN: public RModel{

private:
   
    std::unique_ptr<RFunction> edge_block;
    std::unique_ptr<RFunction> nodes_block;
    std::unique_ptr<RFunction> globals_block;

    std::vector<std::pair<std::string>> edges;
    std::vector<std::string> nodes;
    std::vector<std::string> globals;
    std::vector<std::string> senders;
    std::vector<std::string> receivers;

public:

   //explicit move ctor/assn
   RModel_GNN(RModel_GNN&& other);

   RModel_GNN& operator=(RModel_GNN&& other);

   //disallow copy
   RModel_GNN(const RModel_GNN& other) = delete;
   RModel_GNN& operator=(const RModel_GNN& other) = delete;

   RModel_GNN(){}
   RModel_GNN(std::string name, std::string parsedtime);

   
   ~RModel_GNN(){}}
   ClassDef(RModel_GNN,1);
};

}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RMODEL_GNN