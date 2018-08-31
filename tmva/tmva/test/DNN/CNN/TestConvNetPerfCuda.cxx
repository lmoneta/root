#include <iostream>
#include "TMVA/DNN/Architectures/Cuda.h"

#include "TestConvNetPerf.h"

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

// test here a single operation (e.g. Hadamard)
double testSingleOp1(int nrep, int n0, int n1, int n2) 
{ 
   std::vector<TCudaMatrix<float>> vA; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vA, n0,n1,n2); 
   std::vector<TCudaMatrix<float>> vB; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vB, n0,n1,n2); 

   std::chrono::time_point<std::chrono::system_clock> tstart, tend;
   tstart = std::chrono::system_clock::now();
   for (int i= 0; i < nrep; ++i) {
      for (int j = 0; j < n0; ++j) {
         TCuda<float>::Hadamard(vA[j],vB[j]); 
      }
   }   
   tend = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = tend - tstart;
   return elapsed_seconds.count();       
}      

// test here a single operation (e.g. Hadamard)
double testSingleOp2(int nrep, int n0, int n1, int n2) 
{ 
   std::vector<TCudaMatrix<float>> vA; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vA, n0,n1,n2); 
   std::vector<TCudaMatrix<float>> vB; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vB, n0,n1,n2); 

   for (size_t i = 0; i < vA.size() ; ++i) {
      cudaStream_t s; 
      cudaStreamCreate(&s);       
      vA[i].SetComputeStream(s);
 
   }

   std::chrono::time_point<std::chrono::system_clock> tstart, tend;
   tstart = std::chrono::system_clock::now();
   for (int i= 0; i < nrep; ++i) {
      for (int j = 0; j < n0; ++j) {
         TCuda<float>::Hadamard(vA[j],vB[j]); 
      }
   }   
   //cudaStreamSynchronize(0); 
   cudaDeviceSynchronize(); 

   tend = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = tend - tstart;
   return elapsed_seconds.count();       
}      

double testSingleOp3(int nrep, int n0, int n1, int n2) 
{ 
   std::vector<TCudaMatrix<float>> vA; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vA, 1,n1,n2*n0); 
   std::vector<TCudaMatrix<float>> vB; 
   ConvNetPerfTest<TCuda<float>>::createTensor(vB, 1,n1,n2*n0); 
   TCudaMatrix<float> A0(n1,n2); 

   std::chrono::time_point<std::chrono::system_clock> tstart, tend;
   tstart = std::chrono::system_clock::now();
   for (int i= 0; i < nrep; ++i) {
      //for (int j = 0; j < n0; ++j) {
      TCuda<float>::VHadamard(vA[0],vB[0],A0,n0); 
         //}
   }   
   tend = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = tend - tstart;
   return elapsed_seconds.count();       
}      


void testCNN() { 

   ConvNetPerfTest<TCuda<float>> test;
   
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

void testOP() { 

   std::cout << "Testing Cuda Operation Performances" << std::endl;

   int nrep = 10000; 
   int n0 = 8; 
   int n1 = 64; 
   int n2 = 1; 

   std::cout << "Testing w/0 streams " << std::endl;
   double t1 = testSingleOp1(nrep, n0, n1, n2);
   std::cout <<  "\t---->\tElapsed time = " << t1 << std::endl;

   std::cout << "Testing with streams " << std::endl;
   double t2 = testSingleOp2(nrep, n0, n1, n2);
   std::cout <<  "\t---->\tElapsed time = " << t2 << std::endl;

   std::cout << "Testing with parallel " << std::endl;
   double t3 = testSingleOp3(nrep, n0, n1, n2);
   std::cout <<  "\t---->\tElapsed time = " << t3 << std::endl;

}

int main() { 

   //testCNN(); 
   testOP(); 
}
