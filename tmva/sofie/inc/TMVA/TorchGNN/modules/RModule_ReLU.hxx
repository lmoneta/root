// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * ReLU module.
*/

#ifndef TMVA_SOFIE_RMODULE_RELU_H_
#define TMVA_SOFIE_RMODULE_RELU_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_ReLU: public RModule {
    public:
        /**
         * Construct the ReLU module.
         * 
         * The ReLU operation will be applied element-wise.
         * 
         * @param x The input.
        */
        RModule_ReLU(std::string x) {
            fInputs = {x};
            fArgs = {};
        }

        /** Destruct the module. */
        ~RModule_ReLU() {};

        /**
         * Apply the ReLU operation min(0, x).
        */
        void Forward() {
            const std::vector<float>& x = fInputModules[0] -> GetOutput();
            std::size_t n = x.size();
            fOutput.resize(n);

            for (std::size_t i = 0; i < n; i++) {
                fOutput[i] = (x[i] < 0) ? 0 : x[i];
            }
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is the same as the input shape.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            std::vector<int> shape = fInputModules[0] -> GetShape();
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "ReLU";
        }

        /** 
         * Save parameters.
         * 
         * Does nothing for this module.
         * 
         * @param dir Save directory.
         */
        void SaveParameters([[maybe_unused]] std::string dir) {}

        /**
         * Load saved parameters.
         * 
         * Does nothing for this module.
        */
        void LoadParameters() {}

        /**
         * Load parameters from PyTorch state dictionary.
         * 
         * Does nothing for this module.
         * 
         * @param state_dict The state dictionary.
        */
        void LoadParameters([[maybe_unused]] std::map<std::string, std::vector<float>> state_dict) {}
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_RELU_H_
