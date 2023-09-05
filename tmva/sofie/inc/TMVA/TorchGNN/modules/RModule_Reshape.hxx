// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Reshape module.
*/

#ifndef TMVA_SOFIE_RMODULE_RESHAPE_H_
#define TMVA_SOFIE_RMODULE_RESHAPE_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <algorithm>
#include <stdexcept>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_Reshape: public RModule {
    public:
        /**
         * Construct the reshape module.
         * 
         * Accepts one wildcard dimension (-1).
         * 
         * @param x The input.
         * @param shape The desired shape.
        */
        RModule_Reshape(std::string x, std::vector<int> shape) {
            fOutShape = shape;
            fWildcard = std::find(shape.begin(), shape.end(), -1) - shape.begin();

            // Check shape.
            if (std::any_of(shape.begin(), shape.end(), [](int i){return i == 0;})) {
                throw std::invalid_argument("Dimension cannot be zero.");
            }
            if (std::any_of(shape.begin(), shape.end(), [](int i){return i < -1;})) {
                throw std::invalid_argument("Shape cannot have negative entries (except for the wildcard dimension).");
            }
            if (std::count(shape.begin(), shape.end(), -1) > 1) {
                throw std::invalid_argument("Shape may have at most one wildcard.");
            }

            // Translate shape argument to string.
            std::string shape_arg = "{";
            bool first = true;
            for (int i: shape) {
                if (!first) {
                    shape_arg += ", ";
                } else {
                    first = false;
                }
                shape_arg += i;
            }
            shape_arg += "}";

            fInputs = {x};
            fArgs = {shape_arg};
        }

        /** Destruct the module. */
        ~RModule_Reshape() {};

        /**
         * Simply forward the input.
         * 
         * Reshaping is done through the inferShape() method.
        */
        void Forward() {
            const std::vector<float>& x = fInputModules[0] -> GetOutput();
            fOutput = x;
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is given by the user.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            int cprod = 1;
            for (std::size_t i = 0; i < fOutShape.size(); i++) {
                if (i != fWildcard) {
                    cprod *= fOutShape[i];
                }
            }
            std::vector<int> shape = fOutShape;
            std::vector<float> x = fInputModules[0] -> GetOutput();
            shape[fWildcard] = x.size() / cprod;
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "Reshape";
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
    private:
        std::size_t fWildcard;  // Index of the wildcard dimension.
        std::vector<int> fOutShape;  // Desired output shape.
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_RESHAPE_H_
