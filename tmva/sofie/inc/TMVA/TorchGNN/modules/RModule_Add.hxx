// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Addition module.
*/

#ifndef TMVA_SOFIE_RMODULE_ADD_H_
#define TMVA_SOFIE_RMODULE_ADD_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <gsl/gsl_cblas.h>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_Add: public RModule {
    public:
        /**
         * Construct the addition module.
         * 
         * @param a The first argument.
         * @param b The second argument.
        */
        RModule_Add(std::string a, std::string b) {
            fInputs = {a, b};
            fArgs = {};
        }

        /** Destruct the module. */
        ~RModule_Add() {};

        /**
         * Add the arguments a and b.
        */
        void Forward() {
            const std::vector<float>& a = fInputModules[0] -> GetOutput();
            const std::vector<float>& b = fInputModules[1] -> GetOutput();
            fOutput = b;

            int n = a.size();

            cblas_saxpy(n, 1, a.data(), 1, fOutput.data(), 1);
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is the same as the input shape.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            return fInputModules[0] -> GetShape();
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "Add";
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

#endif  // TMVA_SOFIE_RMODULE_ADD_H_
