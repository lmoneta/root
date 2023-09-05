// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Softmax module.
*/

#ifndef TMVA_SOFIE_RMODULE_SOFTMAX_H_
#define TMVA_SOFIE_RMODULE_SOFTMAX_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <cmath>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_Softmax: public RModule {
    public:
        /**
         * Construct the softmax module.
         * 
         * @param x The input.
        */
        RModule_Softmax(std::string x) {
            fInputs = {x};
            fArgs = {};
        }

        /** Destruct the module. */
        ~RModule_Softmax() {};

        /**
         * Apply the softmax operation exp(x_i) / sum(exp(x_j)).
         * 
         * The sum is taken over the last dimension.
        */
        void Forward() {
            const std::vector<float>& x = fInputModules[0] -> GetOutput();
            int last_dim = fInputModules[0] -> GetShape().back();
            fOutput = x;

            for (std::size_t i = 0; i < x.size(); i += last_dim) {
                float exps[last_dim];
                float sum = 0;
                for (std::size_t j = i; j < i + last_dim; j++) {
                    exps[j - i] = std::exp(x[j]);
                    sum += exps[j - i];
                }
                for (std::size_t j = i; j < i + last_dim; j++) {
                    fOutput[j] = exps[j - i] / sum;
                }
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
            return "Softmax";
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

#endif  // TMVA_SOFIE_RMODULE_SOFTMAX_H_
