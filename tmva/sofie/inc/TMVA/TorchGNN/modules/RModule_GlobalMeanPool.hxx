// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Global mean pooling module.
*/

#ifndef TMVA_SOFIE_RMODULE_GLOBALMEANPOOL_H_
#define TMVA_SOFIE_RMODULE_GLOBALMEANPOOL_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <set>
#include <iostream>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_GlobalMeanPool: public RModule {
    public:
        /**
         * Construct the global mean pooling module.
         * 
         * @param X A node feature matrix of shape (N_1 + ... + N_B,
         * n_features), where N_i denotes the number of nodes from graph i.
         * @param batch The batch vector which assigns each node in x to a
         * specific graph.
        */
        RModule_GlobalMeanPool(std::string X, std::string batch) {
            fInputs = {X, batch};
            fArgs = {};
        }

        /** Destruct the module. */
        ~RModule_GlobalMeanPool() {};

        /**
         * Apply the global mean pooling operation.
        */
        void Forward() {
            const std::vector<float>& x = fInputModules[0] -> GetOutput();
            const std::vector<float>& batch_float = fInputModules[1] -> GetOutput();
            std::vector<int> batch(batch_float.begin(), batch_float.end());

            int n_unique = fOutShape[0];
            int n_features = fOutShape[1];
            
            fOutput = std::vector<float>(n_unique * n_features);

            // TODO: This approach might lead to overflow for large feature
            // values or many nodes. Should we work in logs or use Welford's
            // online algorithm?

            // Sum all entries belonging to same graph.
            std::vector<int> counts(n_unique);
            for (std::size_t i = 0; i < batch.size(); i++) {
                int x_start = i * n_features;
                int out_start = batch[i] * n_features;

                for (int j = 0; j < n_features; j++) {
                    fOutput[out_start + j] += x[x_start + j];
                }
                counts[batch[i]]++;
            }
            
            // Divide all features in a graph by the corresponding number of nodes.
            for (int i = 0; i < n_unique; i++) {
                int out_start = i * n_features;
                for (int j = 0; j < n_features; j++) {
                    fOutput[out_start + j] /= counts[i];
                }
            }
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is (n_unique, n_features), where
         * n_unique denotes the number of graphs in the batch.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            std::vector<float> batch_float = fInputModules[1] -> GetOutput();

            std::vector<int> shape = fInputModules[0] -> GetShape();
            shape[0] = std::set<int>(batch_float.begin(), batch_float.end()).size();
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "GlobalMeanPool";
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

#endif  // TMVA_SOFIE_RMODULE_GLOBALMEANPOOL_H_
