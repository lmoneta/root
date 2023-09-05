/**
 * Tests the TorchGNN ReLU module (RModule_ReLU).
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/ReLUTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_ReLU.hxx"
#include <iostream>

using namespace TMVA::Experimental::SOFIE;

int main() {
    std::vector<float> a = {-2, -1, 0, 1, 2};

    std::vector<float> expected = {0, 0, 0, 1, 2};

    RModel_TorchGNN model = RModel_TorchGNN({"a"}, {{-1}});
    model.AddModule(RModule_ReLU("a"), "out_1");
    std::vector<float> out = model.Forward(a);

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
