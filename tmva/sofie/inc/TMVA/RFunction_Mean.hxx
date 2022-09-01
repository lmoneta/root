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
    
    std::vector<std::string> fInputTensors;
    std::string fOutputTensor;

    void Initialize(std::vector<std::any> InputTensors){
        fInputTensors = std::any_cast<std::vector<GNN_Agg>>(InputTensors);
        function_block->reset(new RModel);
        if(fTarget != FunctionTarget::GLOBALS){
            for(auto& it:InputTensors){
                fInputTensors.emplace_back("Edge_"+it.receiver+"_"+it.sender);
            }
        }
        
        std::unique_ptr<ROperator> op_concat;
        op_concat.reset(new ROperator_Concat<float>(fInputTensors,1,fFuncName+"InputConcat"));
        function_block->AddOperator(std::move(op_concat));

        std::unique_ptr<ROperator> op_reduce_mean;
        op_reduce_mean.reset(new ROperator_Reduce<float,EReduceOpMode::ReduceMean>(0,1,fFuncName+"InputConcat",fOutputTensor));
        function_block->AddOperator(std::move(op_reduce_mean));

        for(int i=0; i<fInputTensors.size(); ++i){
            function_block->AddInputTensorInfo(fInputTensors[i],ETensorType::FLOAT, fGraph->GetTensorShape(fInputTensors[i]));
            function_block->AddInputTensorName(fInputTensors[i]);
        }
        function_block->AddOutputTensorNameList({fOutputTensor});
    }

    std::string Generate(const std::vector<std::any>& InputTensors, int batchSize){
        Initialize(InputTensors);
        function_block->Generate(Options::kGNNComponent, batchSize);
        return function_block->ReturnGenerated();
    }
}

} //SOFIE
} //Experimental
} //TMVA