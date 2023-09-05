/**
 * Tests the saving functionality of TorchGNN.
 * 
 * To run in ROOT terminal:
 * .L path_to_root/tmva/sofie/test/TorchGNN/SaveTest.cxx
 * main()
*/

#include "TMVA/TorchGNN/RModel_TorchGNN.hxx"
#include "TMVA/TorchGNN/modules/RModule_Add.hxx"

using namespace TMVA::Experimental::SOFIE;

int main() {
    std::vector<float> a = {1, 1.5, 2, 2.5};
    std::vector<float> b = {1, 2, 4, 8};
    std::vector<float> c = {0, 0.5, 0, 0.5};

    RModel_TorchGNN model = RModel_TorchGNN({"a", "b", "c"}, {{-1}, {-1}, {-1}});
    model.AddModule(RModule_Add("a", "b"), "out_1");
    model.AddModule(RModule_Add("out_1", "c"), "out_2");
    model.Save("/home/stefan/root-model", "Model", true);

    return 0;
}
