#ifndef TMVA_SOFIE_RFUNCTION_MLP
#define TMVA_SOFIE_RFUNCTION_MLP


#include "TMVA/SOFIE_common.hxx"
#include "TMVA/RFunction.hxx"
#include "TMVA/RModel_GNN.hxx"

#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <limits>
#include <cassert>

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RFunction_Mean: public RFunction{
    
    private:
        std::vector<std::string> fInputTensors;
        std::string fOutputTensor;

    public:
        RFunction_Mean(){}
        void Initialize(std::vector<std::any> InputTensors){
            fInputTensors = std::any_cast<std::vector<GNN_Agg>>(InputTensors);
            function_block->reset(new RModel);
            if(fTarget != FunctionTarget::GLOBALS){
                for(auto& it:InputTensors){
                    fInputTensors.emplace_back("Edge_"+UTILITY::Clean_name(it.receiver)+"_"+UTILITY::Clean_name(it.sender));
                }
            }
            
            std::unique_ptr<ROperator> op_concat;
            op_concat.reset(new ROperator_Concat<float>(fInputTensors,1,fFuncName+"InputConcat"));
            function_block->AddOperator(std::move(op_concat));

            std::unique_ptr<ROperator> op_reduce_mean;
            op_reduce_mean.reset(new ROperator_Reduce<float,EReduceOpMode::ReduceMean>(1,0,fFuncName+"InputConcat",fOutputTensor));
            function_block->AddOperator(std::move(op_reduce_mean));

            for(int i=0; i<fInputTensors.size(); ++i){
                function_block->AddInputTensorInfo(fInputTensors[i],ETensorType::FLOAT, fGraph->GetTensorShape(fInputTensors[i]));
                function_block->AddInputTensorName(fInputTensors[i]);
            }
            function_block->AddOutputTensorNameList({fOutputTensor});
        }

        std::string Generate(const std::string& funcName, const std::vector<std::any>& InputTensors, const std::string& OutputTensor, int batchSize):
        fFuncName(UTILITY::Clean_name(funcName)),fOutputTensor(UTILITY::Clean_name(outputTensor)){
            fOutputTensor = OutputTensor;
            Initialize(InputTensors);
            function_block->Generate(Options::kGNNComponent, batchSize);
            return "\n//--------- GNN_Mean_Agg"+fFuncName+function_block->ReturnGenerated();
        }
}

} //SOFIE
} //Experimental
} //TMVA