#include <iostream>
#include "TMVA/DNN/Architectures/Cpu.h"

#include "TestConvNetPerf.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;



int main()
{

   ConvNetPerfTest<TCpu<float>> test;
   
   std::cout << "Testing CNN Performaces:" << std::endl;
   std::cout << "Input = Ouput  = ( " << test.batchSize << " , " << test.n0 << " , " << test.n1 << "x" << test.n2 << " )\n";
   std::cout << "nrepetitions = " << test.nrep << std::endl << std::endl;
   
   std::cout << "Testing Forward Pass " << std::endl;
   double t1 = test.testConvForwardPass();
   std::cout <<  "\t---->\tElapsed time = " << t1 << "\tTime/batch = " << t1/test.nrep
             << "\t Evts/sec = " << double(test.nrep*test.batchSize)/t1 << std::endl; 

   std::cout << "Testing Backward Pass " << std::endl;
   double t2 = test.testConvBackwardPass();
   std::cout <<  "\t---->\tElapsed time = " << t2 << "\tTime/batch = " << t2/test.nrep
             << "\t Evts/sec = " << double(test.nrep*test.batchSize)/t2 << std::endl; 

}
