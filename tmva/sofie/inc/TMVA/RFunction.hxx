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
        std::unique_ptr<RModel_GNN> fGraph;
        std::unique_ptr<RModel> function_block;
        FunctionTarget fTarget;
        FunctionRelation fRelation;
    public:
        virtual void Initialize(std::vector<std::any> InputTensors) = 0;
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

        void AddInputTensors(std::any){
                if(fType == FunctionType::UPDATE){
                        std::vector<std::size_t> fInputShape = std::any_cast<std::vector<std::size_t>>(inputShape);
                        for(int i=0; i<fInputShape.size(); ++i){
                                function_block->AddInputTensorInfo(fInputTensors[i],ETensorType::FLOAT, fInputShape[i]);
                                function_block->AddInputTensorName(fInputTensors[i]);
                        }
                } else {
                        std::vector<std::vector<std::size_t>> fInputShape = std::any_cast<std::vector<std::vector<std::size_t>>>(inputShape); 
                        for(int i=0; i<fInputShape.size(); ++i){
                                for(int j=0;j<fInputShape[0].size();++j){
                                        function_block->AddInputTensorInfo(fInputTensors[i][j],ETensorType::FLOAT, fInputShape[i][j]);
                                        function_block->AddInputTensorName(fInputTensors[i][j]);
                                }
                        }
                }
        }
        
        void GenerateModel(const std::string& funcName, std::any inputShape){
            fFuncName = UTILITY::Clean_name(funcName);
            Initialize();
            if(inputShape.size() != fInputTensors.size()){
                if(fType == FunctionType::UPDATE)
                        throw std::runtime_error("Passed input shape for GNN Update Function" + fFuncName + "doesn't matches with the input tensor list size");
                else
                        throw std::runtime_error("Passed input shape for GNN Aggregate Function" + fFuncName + "doesn't matches with the input tensor list size");
            }
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
            inferFunc+=").begin());";
            return inferFunc;
        }
};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION