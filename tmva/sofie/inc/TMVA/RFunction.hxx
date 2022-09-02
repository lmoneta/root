#ifndef TMVA_SOFIE_RFUNCTION
#define TMVA_SOFIE_RFUNCTION

#include <any>
#include "TMVA/ROperator.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel;

enum class FunctionType{
        UPDATE=0, AGGREGATE=1
};
enum class FunctionTarget{
        NODES=0, EDGES=1, GLOBALS=2
};
class RFunction: public ROperator{
    std::string fFuncName;
    FunctionType fType;
    FunctionTarget fTarget;
    std::unique_ptr<RModel_GNN> fGraph;
    std::unique_ptr<RModel> function_block;

    virtual void Initialize(std::vector<std::any> InputTensors) = 0;
    virtual std::string Generate(std::string funcName, std::vector<std::any> InputTensors, std::string outputTensor, int batchSize) = 0;
    virtual ~ROperator(){}
};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION