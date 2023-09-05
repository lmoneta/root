/**
 * Tests the TorchGNN global mean pooling module (RModule_GlobalMeanPool).
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/GlobalMeanPoolTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_GlobalMeanPool.hxx"
#include <iostream>

using namespace TMVA::Experimental::SOFIE;

int main() {
    std::vector<float> X = {
        1, 2, 
        -1, -2,
        3, -6,
        1, 2,
        2, 1,
        8, 8
    };
    std::vector<float> batch = {0, 0, 0, 1, 1, 2};

    std::vector<float> expected = {
        1, -2,
        1.5, 1.5,
        8, 8
    };

    RModel_TorchGNN model = RModel_TorchGNN({"X", "batch"}, {{-1, 2}, {-1}});
    model.AddModule(RModule_GlobalMeanPool("X", "batch"), "out_1");
    std::vector<float> out = model.Forward(X, batch);

    std::cout << "Expected:" << std::endl;
    for (float x: expected) {
        std::cout << x << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Actual:" << std::endl;
    for (float x: out) {
        std::cout << x << std::endl;
    }

    return 0;
}
