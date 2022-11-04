#ifndef TMVA_SOFIE_RFUNCTION
#define TMVA_SOFIE_RFUNCTION

#include <any>
#include "TMVA/RModel_GNN.hxx"

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RModel;

enum class GraphType{
        INVALID=0, GNN=1, GraphIndependent=2
};

enum class FunctionType{
        UPDATE=0, AGGREGATE=1
};
enum class FunctionTarget{
        INVALID=0, NODES=1, EDGES=2, GLOBALS=3
};
enum class FunctionReducer{
        INVALID=0, SUM=1, MEAN=2
};
class RFunction{
    protected:
        std::string fFuncName;
        FunctionType fType;
    public:
        RFunction(){}
        virtual ~RFunction(){}
        FunctionType GetFunctionType(){
                return fType;
        }

        RFunction(std::string funcName, FunctionType type):
                fFuncName(UTILITY::Clean_name(funcName)),fType(type){}

};

class RFunction_Update: public RFunction{
        protected:
                std::shared_ptr<RModel> function_block;
                FunctionTarget fTarget;
                GraphType fGraphType;
                std::vector<std::string> fInputTensors;
        public:
        virtual ~RFunction_Update(){}
        RFunction_Update(){}
                RFunction_Update(FunctionTarget target, GraphType gType): fTarget(target), fGraphType(gType){
                        switch(target){
                                case FunctionTarget::EDGES:{
                                        fFuncName = "edge_update";
                                        break;
                                } 
                                case FunctionTarget::NODES: {
                                        fFuncName = "node_update";
                                        break;
                                }
                                case FunctionTarget::GLOBALS: {
                                        fFuncName = "global_update";
                                        break;
                                }
                                default:
                                        throw std::runtime_error("Invalid target for Update function");
                        }
                        fType = FunctionType::UPDATE;
                        function_block = std::make_unique<RModel>(fFuncName);
         
                }
                virtual void AddInitializedTensors(std::any){};
                virtual void Initialize(){};
                void AddInputTensors(std::vector<std::vector<std::size_t>> fInputShape){
                        for(long unsigned int i=0; i<fInputShape.size(); ++i){
                                function_block->AddInputTensorInfo(fInputTensors[i],ETensorType::FLOAT, fInputShape[i]);
                                function_block->AddInputTensorName(fInputTensors[i]);
                        }
                }
                std::shared_ptr<RModel> GetFunctionBlock(){
                        return function_block;
                }

                std::string GenerateModel(std::string filename, long read_pos=0){
                        function_block->SetFilename(filename);
                        function_block->Generate(Options::kGNNComponent,1,read_pos);
                        std::string modelGenerationString;
                        modelGenerationString = "\n//--------- GNN_Update_Function---"+fFuncName+"\n"+function_block->ReturnGenerated();
                        return modelGenerationString;
                }
                std::string Generate(std::vector<std::string> inputPtrs){
                        std::string inferFunc = fFuncName+".infer(";
                        for(auto&it : inputPtrs){
                                inferFunc+=it;
                                inferFunc+=",";
                        }
                        inferFunc.pop_back();
                        inferFunc+=");";
                        return inferFunc;
                }
};

class RFunction_Aggregate: public RFunction{
        protected:
                FunctionReducer fReducer;
        public:
        virtual ~RFunction_Aggregate(){}
        RFunction_Aggregate(){}
                RFunction_Aggregate(FunctionReducer reducer): fReducer(reducer){
                        fType = FunctionType::AGGREGATE;
                }
        virtual std::string GenerateModel() = 0;
        std::string GetFunctionName(){
                return fFuncName;
        }
        FunctionReducer GetFunctionReducer(){
                return fReducer;
        }
        std::string Generate(int num_features, std::vector<std::string> inputTensors){
                std::string inferFunc = fFuncName+"("+std::to_string(num_features)+",{";
                for(auto&it : inputTensors){
                        inferFunc+=it;
                        inferFunc+=",";
                }
                inferFunc.pop_back();
                inferFunc+="});";
                return inferFunc;
        }

};


}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_RFUNCTION
