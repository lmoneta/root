#ifndef TMVA_SOFIE_ROPERATOR_RESHAPE
#define TMVA_SOFIE_ROPERATOR_RESHAPE

#include "TMVA/SOFIE_common.hxx"
#include "TMVA/ROperator.hxx"
#include "TMVA/RModel.hxx"

#include <sstream>

namespace TMVA{
namespace Experimental{
namespace SOFIE{




template <typename T>
class ROperator_Reshape final : public ROperator
{

private:
   std::vector<int_t> fAttrPerm;  

   std::string fNData;        // input data tensor name
   std::string fNShape;       // reshape tensor name
   std::string fNOutput;               // output tensor name
   std::vector<size_t> fShapeData;     // input shape data
   std::vector<size_t> fShapeOutput;   // output shape data

public:

   ROperator_Reshape(){}
   ROperator_Reshape(std::vector<int_t> attr_perm, std::string nameData, std::string nameShape, std::string nameOutput):
      fAttrPerm(attr_perm), fNData(UTILITY::Clean_name(nameData)), fNShape(UTILITY::Clean_name(nameShape)), fNOutput(UTILITY::Clean_name(nameOutput)) {
   }

   ROperator_Reshape(std::string nameData, std::string nameShape, std::string nameOutput)
      : fNData(UTILITY::Clean_name(nameData)), fNShape(UTILITY::Clean_name(nameShape)) , fNOutput(UTILITY::Clean_name(nameOutput))
   {
   }

   std::vector<ETensorType> TypeInference(std::vector<ETensorType> input){
      return input;
   }

   std::vector<std::vector<size_t>> ShapeInference(std::vector<std::vector<size_t>> input){
      if (input.size() != 2) throw std::runtime_error("TMVA SOFIE Reshape Op needs 2 input tensors");
      std::vector<std::vector<size_t>> ret;
      ret.push_back(input[1]);
      return ret;
   }


   void Initialize(RModel& model){
      if (model.CheckIfTensorAlreadyExist(fNData) == false){   //input must be a graph input, or already initialized intermediate tensor
         throw std::runtime_error("TMVA Reshape Op Input Tensor is not found in model");
      }
      fShapeOutput = model.GetTensorShape(fNShape);

      model.AddIntermediateTensor(fNOutput, model.GetTensorType(fNData), fShapeOutput);
   
   }

   std::string Generate(std::string OpName){
      OpName = "op_" + OpName;
      if (fShapeData.empty() || fShapeOutput.empty()){
         throw std::runtime_error("TMVA SOFIE Reshape Op called to Generate without being initialized first");
      }
      // we copy input into output tensor. 
       // output of reshape is same as input
      int length = 1;
      for (auto &i : fShapeOutput) {
          length *= i;
      }
  
      std::string tensor_input = "tensor_" + fNData;
      std::string tensor_output = "tensor_" + fNOutput;
      std::stringstream out;
      out << "\t" << "std::copy(" << tensor_input << "," << tensor_input << " + " << length
          << "," << tensor_output << ")\n";
      
      return out.str(); 

   }


};

}//SOFIE
}//Experimental
}//TMVA


#endif //TMVA_SOFIE_ROPERATOR_RESHAPE
