// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Graph convolution module using Eigen for sparse operations.
 * 
 * This module applies the graph convolution operation X = D^(-1/2) * A *
 * D^(-1/2) * X * Theta.
 * 
 * Its performance is currently worse than the regular GCNConv module, so it is
 * not made available to users.
*/

#ifndef RMODULE_SPARSEGCNCONV_H_
#define RMODULE_SPARSEGCNCONV_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <fstream>
#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCore>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_SparseGCNConv: public RModule {
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
        RModule_SparseGCNConv(std::string x, std::string edge_index, int in_features, int out_features, bool improved=false, bool add_self_loops=true, bool normalize=true, bool bias=true) {
            fInputFeatures = in_features;
            fOutputFeatures = out_features;
            fImprove = improved;
            fSelfLoops = add_self_loops;
            fNormalization = normalize;
            fIncludeBias = bias;
            fEdgeWeights = false;

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
        RModule_SparseGCNConv(std::string x, std::string edge_index, std::string edge_weight, int in_features, int out_features, bool improved=false, bool add_self_loops=true, bool normalize=true, bool bias=true) {
            fInputFeatures = in_features;
            fOutputFeatures = out_features;
            fImprove = improved;
            fSelfLoops = add_self_loops;
            fNormalization = normalize;
            fIncludeBias = bias;
            fEdgeWeights = true;

            if (!fIncludeBias) {
                fB = std::vector<float>(fOutputFeatures);
            }

            fInputs = {x, edge_index, edge_weight};
            fArgs = {std::to_string(in_features), std::to_string(out_features), std::to_string(improved), std::to_string(add_self_loops), std::to_string(normalize), std::to_string(bias)};
        }

        /** Destruct the module. */
        ~RModule_SparseGCNConv() {};

        /**
         * Applies the graph convolution operation to each node.
         * 
         * @returns The updated feature matrix.
        */
        void Forward() {
            const std::vector<float>& X_const = fInputModules[0] -> GetOutput();
            std::vector<float> X = X_const;
            const std::vector<float>& edge_index_f = fInputModules[1] -> GetOutput();

            std::size_t n_nodes = X.size() / fInputFeatures;
            std::size_t n_edges = edge_index_f.size() / 2;
            
            std::vector<float> edge_weight;
            if (fEdgeWeights) {
                const std::vector<float>& edge_weight_f = fInputModules[2] -> GetOutput();
            } else {
                edge_weight = std::vector<float>(n_edges, 1);
            }
            std::vector<float> degree;

            typedef Eigen::Triplet<float> T;
            std::vector<T> edge_list;
            edge_list.reserve(n_edges);
            for(std::size_t i = 0; i < n_edges; i++) {
                int source = edge_index_f[i];
                int target = edge_index_f[i + n_edges];
                edge_list.push_back(T(source, target, edge_weight[i]));
            }
            Eigen::SparseMatrix<float, Eigen::RowMajor> A(n_nodes, n_nodes);
            A.setFromTriplets(edge_list.begin(), edge_list.end());
            
            typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rMatrix;
            if (fNormalization) {
                if (fSelfLoops) {
                    if (fImprove) {
                        A = A + 2 * rMatrix::Identity(n_nodes, n_nodes);
                        degree = std::vector<float>(n_nodes, 2);
                    } else {
                        A = A + rMatrix::Identity(n_nodes, n_nodes);
                        degree = std::vector<float>(n_nodes, 1);
                    }
                } else {
                    degree = std::vector<float>(n_nodes, 0);
                }

                // Loop through edges to get node degrees.
                for (std::size_t i = 0; i < n_edges; i++) {
                    int target = edge_index_f[i + n_edges];
                    degree[target] += edge_weight[i];
                }
            }

            Eigen::Map<Eigen::VectorXf> d(degree.data(), n_nodes);
            d = d.array().sqrt().inverse();
            auto D = d.asDiagonal();
            A = D * A * D;
            
            Eigen::Map<rMatrix> X_m(X.data(), n_nodes, fInputFeatures);
            Eigen::Map<rMatrix> theta(fW.data(), fOutputFeatures, fInputFeatures);
            rMatrix out = A * X_m * theta.transpose();

            if (fIncludeBias) {
                Eigen::Map<Eigen::RowVectorXf> bias(fB.data(), fOutputFeatures);
                out = out.rowwise() + bias;
            }
            fOutput = std::vector<float>(out.data(), out.data() + n_nodes * fOutputFeatures);
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
            return "SparseGCNConv";
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
            std::string del_string = "inc/modules/RModule_SparseGCNConv.hxx";
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
        bool fEdgeWeights;  // True if edge weights are provided.
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

#endif  // RMODULE_SPARSEGCNCONV_H_
