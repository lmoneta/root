// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Input module.
 * 
 * Used internally.
*/

#ifndef TMVA_SOFIE_RMODULE_INPUT_H_
#define TMVA_SOFIE_RMODULE_INPUT_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <algorithm>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_Input: public RModule {
    public:
        /**
         * Construct the input module.
         * 
         * The module stores the input and its shape so that child modules can access it.
         * 
         * @param input_shape The shape of the input. 
        */
        RModule_Input(std::vector<int> input_shape) {
            fInShape = input_shape;
            fWildcard = std::find(input_shape.begin(), input_shape.end(), -1) - input_shape.begin();
            
            fInputs = {};  // No previous inputs to this module.
            fArgs = {};
        }

        /** Destruct the module. */
        ~RModule_Input() {};

        /**
         * Assign input.
         * 
         * @param input The input.
        */
        void SetParams(std::vector<float> input) {
            fOutput = input;
        }

        /**
         * Does nothing for this module.
        */
        void Forward() {}

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is the same as the input shape
         * with an inferred value for the wildcard dimension.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            int cprod = 1;
            for (std::size_t i = 0; i < fInShape.size(); i++) {
                if (i != fWildcard) {
                    cprod *= fInShape[i];
                }
            }
            std::vector<int> shape = fInShape;
            shape[fWildcard] = fOutput.size() / cprod;
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "Input";
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
        std::vector<int> fInShape;  // Input shape.
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_INPUT_H_
