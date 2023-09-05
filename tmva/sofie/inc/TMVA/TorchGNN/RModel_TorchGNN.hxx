// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Header file for PyTorch Geometric models.
 * 
 * Models are created by the user and parameters can then be loaded into each layer.
 * 
 * IMPORTANT: Changes to the format (e.g., namespaces) may affect the emit
 * defined in RModel_TorchGNN.cxx (save).
*/

#ifndef TMVA_SOFIE_RMODEL_TORCHGNN_H_
#define TMVA_SOFIE_RMODEL_TORCHGNN_H_

#include "TMVA/TorchGNN/modules/RModule.hxx"
#include "TMVA/TorchGNN/modules/RModule_Input.hxx"
#include <stdexcept>
#include <iostream>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModel_TorchGNN {
    public:
        /** Model constructor without inputs. */
        RModel_TorchGNN() {}

        /**
         * Model constructor with manual input names.
         * 
         * @param input_names Vector of input names.
         * @param input_shapes Vector of input shapes. Each element may contain
         * at most one wildcard (-1).
        */
        RModel_TorchGNN(std::vector<std::string> input_names, std::vector<std::vector<int>> input_shapes) {
            fInputs = input_names;
            fShapes = input_shapes;

            // Generate input layers.
            for (std::size_t i = 0; i < input_names.size(); i++) {
                // Check shape.
                if (std::any_of(input_shapes[i].begin(), input_shapes[i].end(), [](int j){return j == 0;})) {
                    throw std::invalid_argument("Invalid input shape for input " + input_names[i] + ". Dimension cannot be zero.");
                }
                if (std::any_of(input_shapes[i].begin(), input_shapes[i].end(), [](int j){return j < -1;})) {
                    throw std::invalid_argument("Invalid input shape for input " + input_names[i] + ". Shape cannot have negative entries (except for the wildcard dimension).");
                }
                if (std::count(input_shapes[i].begin(), input_shapes[i].end(), -1) > 1) {
                    throw std::invalid_argument("Invalid input shape for input " + input_names[i] + ". Shape may have at most one wildcard.");
                }
                AddModule(RModule_Input(input_shapes[i]), input_names[i]);
            }
        }

        /**
         * Add a module to the module list.
         * 
         * @param module Module to add.
         * @param name Module name. Defaults to the module type with a count
         * value (e.g., GCNConv_1).
        */
        template<typename T>
        void AddModule(T module, std::string name="") {
            std::string new_name = (name == "") ? std::string(module.GetOperation()) : name;
            if (fModuleCounts[new_name] > 0) {
                // Module exists, so add discriminator and increment count.
                new_name += "_" + std::to_string(fModuleCounts[new_name]);
                fModuleCounts[new_name]++;

                if (name != "") {
                    // Issue warning.
                    std::cout << "WARNING: Module with duplicate name \"" << name << "\" renamed to \"" << new_name << "\"." << std::endl;
                }
            } else {
                // First module of its kind.
                fModuleCounts[new_name] = 1;
            }
            module.SetName(new_name);

            // Initialize the module.
            module.Initialize(fModules, fModuleMap);

            // Add module to the module list.
            fModules.push_back(std::make_shared<T>(module));
            fModuleMap[std::string(module.GetName())] = fModuleCount;
            fModuleCount++;
        }
        
        /**
         * Run the forward function.
         * 
         * @param args Any number of input arguments.
         * @returns The output of the last layer.
        */
        template<class... Types>
        std::vector<float> Forward(Types... args) {
            auto input = std::make_tuple(args...);

            // Instantiate input layers.
            int k = 0;
            std::apply(
                [&](auto&... in) {
                    ((std::dynamic_pointer_cast<RModule_Input>(fModules[k++]) -> SetParams(in)), ...);
                }, input);

            // Loop through and execute modules.
            for (std::shared_ptr<RModule> module: fModules) {
                module -> Execute();
            }

            // Return output of the last layer.
            const std::vector<float>& out_const = fModules.back() -> GetOutput();
            std::vector<float> out = out_const;
            return out;
        }

        /**
         * Load parameters from PyTorch state dictionary for all modules.
         * 
         * @param state_dict The state dictionary.
        */
        void LoadParameters(std::map<std::string, std::vector<float>> state_dict) {
            for (std::shared_ptr<RModule> module: fModules) {
                module -> LoadParameters(state_dict);
            }
        }

        /**
         * Load saved parameters for all modules.
        */
        void LoadParameters() {
            for (std::shared_ptr<RModule> module: fModules) {
                module -> LoadParameters();
            }
        }

        /**
         * Save the model as standalone inference code.
         * 
         * @param path Path to save location.
         * @param name Model name.
         * @param overwrite True if any existing directory should be
         * overwritten. Defaults to false.
        */
        void Save(std::string path, std::string name, bool overwrite=false);
    private:
        /**
         * Get a timestamp.
         * 
         * @returns The timestamp in string format.
        */
        static std::string GetTimestamp() {
            time_t rawtime;
            struct tm * timeinfo;
            char timestamp [80];
            time(&rawtime);
            timeinfo = localtime(&rawtime);
            strftime(timestamp, 80, "Timestamp: %d-%m-%Y %T.", timeinfo);
            return timestamp;
        }

        /**
         * Write the methods to create a self-contained package.
         * 
         * @param dir Directory to save to.
         * @param name Model name.
         * @param timestamp Timestamp.
        */
        void WriteMethods(std::string dir, std::string name, std::string timestamp);

        /**
         * Write the model to a file.
         * 
         * @param dir Directory to save to.
         * @param name Model name.
         * @param timestamp Timestamp.
        */
        void WriteModel(std::string dir, std::string name, std::string timestamp);

        /**
         * Write the CMakeLists file.
         * 
         * @param dir Directory to save to.
         * @param name Model name.
         * @param timestamp Timestamp.
        */
        void WriteCMakeLists(std::string dir, std::string name, std::string timestamp);

        std::vector<std::string> fInputs;  // Vector of input names.
        std::vector<std::vector<int>> fShapes;  // Vector of input shapes.
        std::map<std::string, int> fModuleCounts;  // Map from module name to number of occurrences.
        std::vector<std::shared_ptr<RModule>> fModules;  // Vector containing the modules.
        std::map<std::string, int> fModuleMap;  // Map from module name to module index (in modules).
        int fModuleCount = 0;  // Number of modules.
};

}  // SOFIE.
}  // Experimental.
}  // TMVA.

#endif  // TMVA_SOFIE_RMODEL_TORCHGNN_H_
