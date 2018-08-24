
using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

#include "TestConvNet.h"

#include <chrono>

size_t batchSize = 32; 
int nrep = 1000; 

// conv layer output dimension matrix 
size_t nout1 = 12;   
size_t nout2 = 32*32;   

template<typename Architecture>
void createTensor(std::vector<size_t n , size_t m)
{

   using Matrix_t = typename Architecture::Matrix_t;


template<typename Architecture>
bool testEvalFunction()
{

   using Matrix_t = typename Architecture::Matrix_t;

   std::vector<Matrix_t> & vec; 
   for (size_t i = 0; i < batchSize; ++i)  {
      vec.emplace_back(nout1,nout2); 
      Architecture::InitializeGlorotNormal( vec[i] ); 
   }
   
   for (int i = 0; i < 

}
