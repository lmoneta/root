// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Graph convolution module.
 * 
 * This module applies the graph convolution operation X = D^(-1/2) * A *
 * D^(-1/2) * X * Theta.
*/

#ifndef TMVA_SOFIE_RMODULE_GCNCONV_H_
#define TMVA_SOFIE_RMODULE_GCNCONV_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <gsl/gsl_cblas.h>
#include <cmath>
#include <fstream>
#include <iostream>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_GCNConv: public RModule {
    public:
        /**
         * Construct the graph convolution module without edge weights.
         * 
         * @param x The input feature matrix of shape (n_nodes, in_features).
         * @param edge_index The edge indices matrix of shape (2, n_edges).
         * @param in_features The size of each input sample.
         * @param out_features The size of each output sample.
         * @param improved True if self-loops should have a weight of two (A = A +
         * 2I). Defaults to false.
         * @param add_self_loops True if self-loops should be added (A = A + I).
         * Defaults to true.
         * @param normalize True if self-loops should be added and symmetric
         * normalization should be computed on the fly.
         * @param bias True if a bias is included. Defaults to true.
        */
        RModule_GCNConv(std::string x, std::string edge_index, int in_features, int out_features, bool improved=false, bool add_self_loops=true, bool normalize=true, bool bias=true) {
            fInputFeatures = in_features;
            fOutputFeatures = out_features;
            fImprove = improved;
            fSelfLoops = add_self_loops;
            fNormalization = normalize;
            fIncludeBias = bias;
            fUseEdgeWeights = false;

            if (!fIncludeBias) {
                fB = std::vector<float>(fOutputFeatures);
            }

            fInputs = {x, edge_index};
            fArgs = {std::to_string(in_features), std::to_string(out_features), std::to_string(improved), std::to_string(add_self_loops), std::to_string(normalize), std::to_string(bias)};
        }

        /**
         * Construct the graph convolution module with edge weights.
         * 
         * @param x The input feature matrix of shape (n_nodes, in_features).
         * @param edge_index The edge indices matrix of shape (2, n_edges).
         * @param edge_weight The edge weights vector of shape (n_edges).
         * @param in_features The size of each input sample.
         * @param out_features The size of each output sample.
         * @param improved True if self-loops should have a weight of two (A = A +
         * 2I). Defaults to false.
         * @param add_self_loops True if self-loops should be added (A = A + I).
         * Defaults to true.
         * @param normalize True if self-loops should be added and symmetric
         * normalization should be computed on the fly.
         * @param bias True if a bias is included. Defaults to true.
        */
        RModule_GCNConv(std::string x, std::string edge_index, std::string edge_weight, int in_features, int out_features, bool improved=false, bool add_self_loops=true, bool normalize=true, bool bias=true) {
            fInputFeatures = in_features;
            fOutputFeatures = out_features;
            fImprove = improved;
            fSelfLoops = add_self_loops;
            fNormalization = normalize;
            fIncludeBias = bias;
            fUseEdgeWeights = true;

            if (!fIncludeBias) {
                fB = std::vector<float>(fOutputFeatures);
            }

            fInputs = {x, edge_index, edge_weight};
            fArgs = {std::to_string(in_features), std::to_string(out_features), std::to_string(improved), std::to_string(add_self_loops), std::to_string(normalize), std::to_string(bias)};
        }

        /** Destruct the module. */
        ~RModule_GCNConv() {};

        /**
         * Applies the graph convolution operation to each node.
        */
        void Forward() {
            const std::vector<float>& X = fInputModules[0] -> GetOutput();
            const std::vector<float>& edge_index = fInputModules[1] -> GetOutput();

            std::size_t n_nodes = X.size() / fInputFeatures;
            std::size_t n_edges = edge_index.size() / 2;
            
            std::vector<float> edge_weight;
            if (fUseEdgeWeights) {
                edge_weight = fInputModules[2] -> GetOutput();
            } else {
                edge_weight = std::vector<float>(n_edges, 1);
            }
            std::vector<float> X_agg;
            std::vector<float> degree;
            
            if (fNormalization) {
                if (fSelfLoops) {
                    if (fImprove) {
                        degree = std::vector<float>(n_nodes, 2);
                    } else {
                        degree = std::vector<float>(n_nodes, 1);
                    }
                } else {
                    degree = std::vector<float>(n_nodes, 0);
                }

                // Loop through edges to get node degrees.
                for (std::size_t i = 0; i < n_edges; i++) {
                    int target = edge_index[i + n_edges];
                    degree[target] += edge_weight[i];
                }
            }

            if (fNormalization && fSelfLoops) {
                // Add self loops.
                X_agg = X;
                int self_weight;
                if (fImprove) {
                    self_weight = 2;
                } else {
                    self_weight = 1;
                }
                for (std::size_t i = 0; i < n_nodes; i++) {
                    for (int j = 0; j < fInputFeatures; j++) {
                        X_agg[i * fInputFeatures + j] *= self_weight / degree[i];
                    }
                }
            } else {
                // Set X_agg to zero.
                X_agg = std::vector<float>(X.size());
            }

            // Loop through the edges to aggregate information from neighboring
            // nodes.
            for (std::size_t i = 0; i < n_edges; i++) {
                int source = edge_index[i];
                int target = edge_index[i + n_edges];
                
                int x_start = source * fInputFeatures;
                int x_agg_start = target * fInputFeatures;
                float norm = edge_weight[i] / std::sqrt(degree[source] * degree[target]);
                for (int j = 0; j < fInputFeatures; j++) {
                    if (fNormalization) {
                        X_agg[x_agg_start + j] += norm * X[x_start + j];
                    } else {
                        X_agg[x_agg_start + j] += edge_weight[i] * X[x_start + j];
                    }
                }
            }

            // Update node features.
                // Perform matrix multiplications (Y = X_agg * W^T + 1b^T :=
                // X_agg W^T + B).
                    // X_agg, shape (m, k) -> (n_nodes, input_features).
                    // W, shape (n, k) -> (output_features, input_features).
                    // B, shape (m, n) -> (n_nodes, output_features).
            int m = n_nodes;
            int k = fInputFeatures;
            int n = fOutputFeatures;

            fOutput = std::vector<float>(n_nodes * fOutputFeatures, 0);  // cblas sets (B = X_agg * W^T + B).
            if (fIncludeBias) {
                // Construct bias matrix (B = 1b^T).
                std::vector<float> one(n_nodes, 1);
                cblas_sger(CblasRowMajor, n_nodes, fOutputFeatures, 1, one.data(), 1, fB.data(), 1, fOutput.data(), fOutputFeatures);
            }

            // Perform matrix multiplication (Y = X_agg * W^T + B).
            cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1, X_agg.data(), k, fW.data(), k, 1, fOutput.data(), n);
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is the same as the shape of the
         * input feature matrix but with out_features on the last dimension 
         * instead of in_features.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            std::vector<int> shape = fInputModules[0] -> GetShape();
            shape.back() = fOutputFeatures;
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "GCNConv";
        }

        /**
         * Set the weights.
         * 
         * @param weights The weight matrix.
        */
        void SetWeights(std::vector<float> weights) {fW = weights;}

        /**
         * Set the biases.
         * 
         * @param biases The bias vector.
        */
        void SetBiases(std::vector<float> biases) {fB = biases;}

        /** 
         * Save parameters.
         * 
         * @param dir Save directory.
         */
        void SaveParameters(std::string dir) {
            // Save weights.
            std::string fdir = dir + "/" + fName + "_lin_weight.dat";
            std::ofstream outfile = std::ofstream(fdir, std::ios::out | std::ios::binary);
            outfile.write(reinterpret_cast<char*>(&fW[0]), fW.size() * sizeof(float));
            
            if (fIncludeBias) {
                // Save biases.
                fdir = dir + "/" + fName + "_bias.dat";
                outfile = std::ofstream(fdir, std::ios::out | std::ios::binary);
                outfile.write(reinterpret_cast<char*>(&fB[0]), fB.size() * sizeof(float));
            }
            outfile.close();
        }

        /**
         * Load saved parameters.
        */
        void LoadParameters() {
            std::string dir = __FILE__;
            std::string del_string = "inc/modules/RModule_GCNConv.hxx";
            dir.replace(dir.find(del_string), del_string.size(), "params/");

            // Load weights.
            std::string param_dir = dir + fName + "_lin_weight.dat";
            std::ifstream infile = std::ifstream(param_dir, std::ios::in | std::ios::binary);
            fW = std::vector<float>(fInputFeatures * fOutputFeatures);
            infile.read(reinterpret_cast<char*>(&fW[0]), fW.size() * sizeof(float));
            infile.close();

            if (fIncludeBias) {
                // Load biases.
                param_dir = dir + fName + "_bias.dat";
                infile = std::ifstream(param_dir, std::ios::in | std::ios::binary);
                fB = std::vector<float>(fOutputFeatures);
                infile.read(reinterpret_cast<char*>(&fB[0]), fB.size() * sizeof(float));
                infile.close();
            }
        }

        /**
         * Load parameters from PyTorch state dictionary.
         * 
         * @param state_dict The state dictionary.
        */
        void LoadParameters(std::map<std::string, std::vector<float>> state_dict) {
            if (auto search = state_dict.find(fName + ".lin.weight"); search != state_dict.end()) {
                fW = state_dict[fName + ".lin.weight"];
            } else {
                std::cout << "WARNING: Weights for module " << fName << " not found." << std::endl;
            }
            
            if (fIncludeBias) {
                if (auto search = state_dict.find(fName + ".bias"); search != state_dict.end()) {
                    fB = state_dict[fName + ".bias"];
                } else {
                    std::cout << "WARNING: Biases for module " << fName << " not found." << std::endl;
                }
            }
        }
    private:
        int fInputFeatures;  // The size of each input sample.
        int fOutputFeatures;  // The size of each output sample.
        bool fUseEdgeWeights;  // True if edge weights are provided.
        bool fImprove;  // True if self-loops should get a weight of two.
        bool fSelfLoops;  // True if self-loops should be added.
        bool fNormalization;  // True if edge weights should be normalized.
        bool fIncludeBias;  // True if a bias is included.
        std::vector<float> fW;  // Weight matrix W.
        std::vector<float> fB;  // Bias vector b.
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_GCNCONV_H_
