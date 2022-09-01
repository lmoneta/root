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
#include <any>

namespace TMVA{
namespace Experimental{
namespace SOFIE{

class RFunction_MLP: public RFunction{
    Int_t fNumLayers;          // Number of Layers in MLP
    bool useActivation;       // if True, ReLU is used as activation for every layer of the MLP

    std::vector<std::string> fInputTensors;
    std::vector<std::string> fKernelTensors;
    std::vector<std::string> fBiasTensors;
    std::string fOutputTensor;

    void Initialize(std::vector<std::any> InputTensors){
        
        function_block->reset(new RModel);
        
        if(fTarget == FunctionTarget::EDGES){
            fInputTensors.emplace_back("Edge_"+std::any_cast<std::string>(InputTensors[0])+"_"+std::any_cast<std::string>(InputTensors[1]))
        } else if(fTarget == FunctionTarget::NODES){
            auto nodes_agg_data = std::any_cast<std::vector<GNN::GNN_Agg>>(InputTensors[0]);
            for(auto& it:nodes_agg_data){
                fInputTensors.emplace_back("Edge_"+it.receiver+"_"+it.sender);
            }
        }

        for(int i=1; i<InputTensors.size();++i){
            fInputTensors.emplace_back(std::any_cast<std::string>(InputTensors[i]));
        }

        std::unique_ptr<ROperator> op_concat;
        op_concat.reset(new ROperator_Concat<float>(fInputTensors,0,fFuncName+"InputConcat"));
        function_block->AddOperator(std::move(op_concat));
        
        std::unique_ptr<ROperator> op_gemm;
        std::string fGemmInput = fFuncName+"InputConcat";
        for(int i=0; i<fNumLayers; ++i){
            op_gemm.reset(new ROperator_Gemm<float>(1.0,1.0,0,0,fGemmInput,fKernelTensors[i],fBiasTensors[i],fFuncName+"Gemm"+i));
            function_block->AddOperator(std::move(op_gemm));
            fGemmInput = fFuncName+"Gemm"+i;
            if(useActivation){
                std::unique_ptr<ROperator> op_relu;
                op_relu.reset(new ROperator_Relu<float>(fFuncName+"Gemm"+i, fFuncName+"Relu"+i));
                function_block->AddOperator(std::move(op_relu));
                fGemmInput = fFuncName+"Relu"+i;
            }
        }
        function_block->AddBlasRoutines({"Gemm", "Gemv"});  // for Gemm operation

        // assuming all the linear layers has a kernel and a bias initialized tensors
        for(int i=0;i<fKernelTensors.size();++i){
            function_block->AddInitializedTensor(fKernelTensors[i],ETensorType::FLOAT,fGraph->GetTensorShape(fKernelTensors[i]),fGraph->GetInitializedTensorData(fKernalTensors[i]));
            function_block->AddInitializedTensor(fBiasTensors[i],ETensorType::FLOAT,fGraph->GetTensorShape(fBiasTensors[i]),fGraph->GetInitializedTensorData(fBiasTensors[i]));
        }
        
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

} // SOFIE
} // Experimental
} // TMVA

#endif //TMVA_SOFIE_RFUNCTION_MLP