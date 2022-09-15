#ifndef TMVA_SOFIE_RFUNCTION
#define TMVA_SOFIE_RFUNCTION

#include <any>
#include "TMVA/RModel_GNN.hxx"
#include "TMVA/ROperator.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel;
class RModel_GNN;

enum class FunctionType{
        UPDATE=0, AGGREGATE=1
};
enum class FunctionTarget{
        INVALID=0, NODES=1, EDGES=2, GLOBALS=3
};
enum class FunctionRelation{
        INVALID=0, NODES_GLOBALS=1, EDGES_GLOBALS=2, EDGES_NODES=3
};
class RFunction: public ROperator{
    protected:
        std::string fFuncName;
        FunctionType fType;
        std::unique_ptr<RModel> function_block;
        FunctionTarget fTarget;
        FunctionRelation fRelation;
    public:
        virtual void Initialize() = 0;
        virtual ~RFunction(){}
        FunctionType GetFunctionType(){
                return fType;
        }
        FunctionTarget GetFunctionTarget(){
                return fTarget;
        }
        FunctionRelation GetFunctionRelation(){
                return fRelation;
        }
        std::unique_ptr<RModel> GetFunctionBlock(){
                return std::move(function_block);
        }

        RFunction(FunctionType Type,FunctionTarget target, FunctionRelation relation):
                fType(Type), fTarget(target), fRelation(relation){}

        virtual void AddInputTensors(std::any inputShape) = 0;
        
        std::string GenerateModel(const std::string& funcName, std::any inputShape){
            fFuncName = UTILITY::Clean_name(funcName);
            Initialize();
            AddInputTensors(inputShape);
            function_block->Generate(Options::kGNNComponent);
            std::string modelGenerationString;
            if(fType == FunctionType::UPDATE)
                modelGenerationString = "\n//--------- GNN_Update_Function"+fFuncName+"\n"+function_block->ReturnGenerated();
            else        
                modelGenerationString = "\n//--------- GNN_Aggregate_Function"+fFuncName+"\n"+function_block->ReturnGenerated();
            return modelGenerationString;
        }

        std::string Generate(std::vector<std::string> inputPtrs){
            std::string inferFunc = fFuncName+"::infer(";
            for(auto&it : inputPtrs){
                inferFunc+=it;
                inferFunc+=",";
            }
            inferFunc+=");";
            return inferFunc;
        }
};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION