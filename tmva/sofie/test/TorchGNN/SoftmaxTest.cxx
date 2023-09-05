/**
 * Tests the TorchGNN softmax module (RModule_Softmax).
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/SoftmaxTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_Softmax.hxx"
#include <iostream>

using namespace TMVA::Experimental::SOFIE;

int main() {
    // 4x3.
    std::vector<float> a = {
        0, 3, -2,
        0, 3, -2,
        1, 2, 3,
        10, -1, -1
    };

    // 4x3.
    std::vector<float> expected = {
        0.0471, 0.9465, 0.0064,
        0.0471, 0.9465, 0.0064,
        0.0900, 0.2447, 0.6652,
        1.0000, 0.0000, 0.0000
    };

    RModel_TorchGNN model = RModel_TorchGNN({"a"}, {{-1, 3}});
    model.AddModule(RModule_Softmax("a"), "out_1");
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
