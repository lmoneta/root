// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Linear module.
 * 
 * This module applies a linear transformation (Ax + b) to the data. For
 * matrix inputs of shape (n_obs, n_features), it applies the transformation 
 * (XA^T + 1b^T).
*/

#ifndef TMVA_SOFIE_RMODULE_LINEAR_H_
#define TMVA_SOFIE_RMODULE_LINEAR_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include <gsl/gsl_cblas.h>
#include <fstream>
#include <iostream>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule_Linear: public RModule {
    public:
        /**
         * Construct the linear module.
         * 
         * @param x The input.
         * @param in_features The size of each input sample.
         * @param out_features The size of each output sample.
         * @param bias True if a bias is included. Defaults to true.
        */
        RModule_Linear(std::string x, int in_features, int out_features, bool bias=true) {
            fInputFeatures = in_features;
            fOutputFeatures = out_features;
            fIncludeBias = bias;

            if (!fIncludeBias) {
                fB = std::vector<float>(fOutputFeatures);
            }

            fInputs = {x};
            fArgs = {std::to_string(in_features), std::to_string(out_features), std::to_string(bias)};
        }

        /** Destruct the module. */
        ~RModule_Linear() {};

        /**
         * Applies the linear transformation (y = Ax + b) to each element in the
         * input.
        */
        void Forward() {
            const std::vector<float>& in = fInputModules[0] -> GetOutput();
            

            if (fRowDim > 1) {
                // Perform matrix multiplications (Y = XA^T + 1b^T := XA^T + B).
                // X, shape (m, k) -> (row_dim, input_features).
                // A, shape (n, k) -> (output_features, input_features).
                // B, shape (m, n) -> (row_dim, output_features).
                int m = fRowDim;
                int k = fInputFeatures;
                int n = fOutputFeatures;
                for (std::size_t i = 0; i < in.size(); i += fRowDim * fInputFeatures) {
                    // Get input matrix.
                    std::vector<float> X(in.begin() + i, in.begin() + i + fRowDim * fInputFeatures);

                    fOutput = std::vector<float>(fRowDim * fOutputFeatures, 0);  // cblas sets (B = XA^T + B).
                    if (fIncludeBias) {
                        // Construct bias matrix (B = 1b^T).
                        std::vector<float> one(fRowDim, 1);
                        cblas_sger(CblasRowMajor, fRowDim, fOutputFeatures, 1, one.data(), 1, fB.data(), 1, fOutput.data(), fOutputFeatures);
                    }

                    // Perform matrix multiplication (Y = XA^T + B).
                    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1, X.data(), k, fA.data(), k, 1, fOutput.data(), n);
                }
            } else {
                fOutput.clear();
                fOutput.reserve(fNumOut);

                // Perform matrix-vector multiplications (y = Ax + b).
                // A, shape (n, m) -> (output_features, input_features).
                // x, shape (m, 1) -> (input_features, 1).
                // b, shape (n, 1) -> (output_features, 1).
                int m = fInputFeatures;
                int n = fOutputFeatures;
                for (std::size_t i = 0; i < in.size(); i += fInputFeatures) {
                    // Get input vector.
                    std::vector<float> x(in.begin() + i, in.begin() + i + fInputFeatures);

                    // Copy bias.
                    std::vector<float> y(fB);

                    // Perform matrix-vector multiplication (y = Ax + b).
                    cblas_sgemv(CblasRowMajor, CblasNoTrans, n, m, 1, fA.data(), m, x.data(), 1, 1, y.data(), 1);
                    for (float elem: y) {  // cblas sets (b = Ax + b).
                        fOutput.push_back(elem);
                    }
                }
            }
        }

        /**
         * Infer the output shape.
         * 
         * For this module, the output shape is the same as the input shape but
         * with out_features on the last dimension instead of in_features.
         * 
         * @returns The output shape.
        */
        std::vector<int> InferShape() {
            std::vector<int> shape = fInputModules[0] -> GetShape();
            shape.back() = fOutputFeatures;

            if (shape.size() > 1) {
                fRowDim = shape[shape.size() - 2];
            }

            fNumOut = 1;
            for (int dim: shape) {
                fNumOut *= dim;
            }
            return shape;
        }

        /**
         * Get the operation.
         * 
         * @returns The name of the operation.
        */
        std::string_view GetOperation() {
            return "Linear";
        }

        /**
         * Set the weights.
         * 
         * @param weights The weight matrix.
        */
        void SetWeights(std::vector<float> weights) {fA = weights;}

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
            std::string fdir = dir + "/" + fName + "_weight.dat";
            std::ofstream outfile = std::ofstream(fdir, std::ios::out | std::ios::binary);
            outfile.write(reinterpret_cast<char*>(&fA[0]), fA.size() * sizeof(float));
            
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
            std::string del_string = "inc/modules/RModule_Linear.hxx";
            dir.replace(dir.find(del_string), del_string.size(), "params/");

            // Load weights.
            std::string param_dir = dir + fName + "_weight.dat";
            std::ifstream infile = std::ifstream(param_dir, std::ios::in | std::ios::binary);
            fA = std::vector<float>(fInputFeatures * fOutputFeatures);
            infile.read(reinterpret_cast<char*>(&fA[0]), fA.size() * sizeof(float));
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
            if (auto search = state_dict.find(fName + ".weight"); search != state_dict.end()) {
                fA = state_dict[fName + ".weight"];
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
        bool fIncludeBias;  // True if a bias is included.
        int fRowDim = 1;  // Size of the second to last input dimension.
        int fNumOut;  // The number of output elements.
        std::vector<float> fA;  // Weight matrix A.
        std::vector<float> fB;  // Bias vector b.
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_LINEAR_H_
