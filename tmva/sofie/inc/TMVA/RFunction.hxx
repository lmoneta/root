#ifndef TMVA_SOFIE_RFUNCTION
#define TMVA_SOFIE_RFUNCTION

#include <any>
#include "TMVA/RModel_GNN.hxx"

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
class RFunction{
    protected:
        std::string fFuncName;
        FunctionType fType;
        std::unique_ptr<RModel> function_block;
    public:
        virtual void Initialize() = 0;
        virtual ~RFunction(){}
        
        FunctionType GetFunctionType(){
                return fType;
        }
        std::unique_ptr<RModel> GetFunctionBlock(){
                return std::move(function_block);
        }

        RFunction(std::string funcName, FunctionType type):
                fFuncName( UTILITY::Clean_name(funcName)),fType(type){
                        function_block.reset(new RModel(fFuncName));   
        }

        virtual void AddInputTensors(std::any inputShape) = 0;

        virtual void AddInitializedTensor(std::any);
        
        std::string GenerateModel(std::any inputShape){
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

class RFunction_Update: public RFunction{
        protected:
                FunctionTarget fTarget;
                std::vector<std::string> fInputTensors;
        public:
                RFunction_Update(std::string funcName, FunctionTarget target):RFunction(funcName,FunctionType::UPDATE), fTarget(target){}

                void AddInputTensors(std::any inputShape){
                        std::vector<std::vector<std::size_t>> fInputShape = std::any_cast<std::vector<std::vector<std::size_t>>>(inputShape);
                        for(long unsigned int i=0; i<fInputShape.size(); ++i){
                                function_block->AddInputTensorInfo(fInputTensors[i],ETensorType::FLOAT, fInputShape[i]);
                                function_block->AddInputTensorName(fInputTensors[i]);
                        }
                }
};

class RFunction_Aggregate: public RFunction{
        protected:
                FunctionRelation fRelation;
                std::vector<std::vector<std::string>> fInputTensors;
        public:
                RFunction_Aggregate(std::string funcName, FunctionRelation relation):RFunction(funcName,FunctionType::AGGREGATE), fRelation(relation){}
                void AddInputTensors(std::any inputShape){
                        std::vector<std::vector<std::vector<std::size_t>>> fInputShape = std::any_cast<std::vector<std::vector<std::vector<std::size_t>>>>(inputShape); 
                                for(long unsigned int i=0; i<fInputShape.size(); ++i){
                                        for(long unsigned int j=0;j<fInputShape[0].size();++j){
                                                function_block->AddInputTensorInfo(fInputTensors[i][j],ETensorType::FLOAT, fInputShape[i][j]);
                                                function_block->AddInputTensorName(fInputTensors[i][j]);
                                        }
                                }
                }
};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION