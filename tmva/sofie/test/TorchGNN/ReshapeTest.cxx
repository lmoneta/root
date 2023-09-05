/**
 * Tests the TorchGNN reshape module (RModule_Reshape).
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/ReshapeTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_Cat.hxx"
#include "TMVA/TorchGNN/modules/RModule_Reshape.hxx"
#include <iostream>

using namespace TMVA::Experimental::SOFIE;

int main() {
    // 3x2.
    std::vector<float> a = {
        1, 2,
        3, 4,
        5, 6
    };
    // 1x3.
    std::vector<float> b = {
        -1, -2, -3
    };

    std::vector<float> expected = {
        1, 2, -1,
        3, 4, -2,
        5, 6, -3
    };

    RModel_TorchGNN model = RModel_TorchGNN({"a", "b"}, {{-1, 2}, {-1, 3}});
    std::vector<int> shape = {-1, 1};
    model.AddModule(RModule_Reshape("b", shape), "out_1");
    model.AddModule(RModule_Cat("a", "out_1", 1), "out_2");
    std::vector<float> out = model.Forward(a, b);

    std::cout << "Expected:" << std::endl;
    for (float x: expected) {
        std::cout << x << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Actual:" << std::endl;
    for (float x: out) {
        std::cout << x << std::endl;
    }
    std::cout << std::endl;

    model = RModel_TorchGNN({"a", "b"}, {{-1, 2}, {-1, 3}});
    model.AddModule(RModule_Cat("a", "b", 1), "out_1");

    std::cout << "Invalid shapes should throw an error." << std::endl;
    std::cout << std::endl;

    try {
        std::vector<float> out = model.Forward(a, b);
    } catch (std::string error) {
        std::cout << error << std::endl;;
    }

    return 0;
}
