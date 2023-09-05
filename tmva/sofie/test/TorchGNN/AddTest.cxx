/**
 * Tests the TorchGNN addition module (RModule_Add).
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/AddTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_Add.hxx"
#include <iostream>

using namespace TMVA::Experimental::SOFIE;

int main() {
    std::vector<float> a = {1, 1.5, 2, 2.5};
    std::vector<float> b = {1, 2, 4, 8};
    std::vector<float> c = {0, 0.5, 0, 0.5};

    std::vector<float> expected = {2, 3.5, 6, 10.5};

    RModel_TorchGNN model = RModel_TorchGNN({"a", "b"}, {{-1}, {-1}});
    model.AddModule(RModule_Add("a", "b"), "out_1");
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

    expected = {2, 4, 6, 11};

    model = RModel_TorchGNN({"a", "b", "c"}, {{-1}, {-1}, {-1}});
    model.AddModule(RModule_Add("a", "b"), "out_1");
    model.AddModule(RModule_Add("out_1", "c"), "out_2");
    out = model.Forward(a, b, c);

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

    expected = {2, 4, 6, 11};

    model = RModel_TorchGNN({"a", "a", "c"}, {{-1}, {-1}, {-1}});
    model.AddModule(RModule_Add("a", "a_1"), "out_1");
    model.AddModule(RModule_Add("out_1", "c"), "out_2");
    out = model.Forward(a, b, c);

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
