// @(#)root/tmva/sofie:$Id$
// Author: Stefan van Berkum

/**
 * Base class for RModule objects used in GNNs.
 * 
 * Modules define the operations that can be performed in a forward pass. They
 * can be layers, activations, or generic operations. Each module defines a
 * forward method that updates fOutput (the output of the last call).
 * 
 * IMPORTANT: Besides the virtual methods, each RModule should assign its inputs
 * to the class variable "fInputs" and other arguments to the class variable
 * "fArgs" (in string format).
 * IMPORTANT: Changes to the format (e.g., namespaces) may affect the emit
 * defined in RModel_TorchGNN.cxx (save). To be safe, new modules should closely
 * follow the format of the exisiting modules (i.e., copy-paste and edit).
*/

#ifndef TMVA_SOFIE_RMODULE_H_
#define TMVA_SOFIE_RMODULE_H_

#include <vector>
#include <memory>
#include <map>
#include <string>

namespace TMVA {
namespace Experimental {
namespace SOFIE {

class RModule {
    public:
        virtual ~RModule() {};

        /**
         * Initialize the RModule by binding it to its input modules.
         * 
         * @param module_list Vector containing the modules.
         * @param module_map Map from module name to module index (in modules).
        */
        void Initialize(std::vector<std::shared_ptr<RModule>> module_list, std::map<std::string, int> module_map) {
            for (std::string input: fInputs) {
                fInputModules.push_back(module_list[module_map[input]]);
            }
        }

        /**
         * Execute the module.
         * 
         * This triggers the module's forward method and stores the output and
         * its shape.
        */
        void Execute() {
            fOutShape = InferShape();  // Infer shape on the fly. TODO: Test impact on performance and possibly execute once for static models.
            Forward();
        }

        /**
         * Change this module's name.
         * 
         * @param new_name New module name.
        */
        void SetName(std::string new_name) {fName = new_name;}

        /**
         * Get this module's name.
         * 
         * @returns The module name.
        */
        std::string_view GetName() {return fName;}

        /**
         * Get the module inputs.
         * 
         * @returns The inputs for this module.
        */
        std::vector<std::string> GetInputs() {return fInputs;}

        /**
         * Get the module arguments.
        */
        std::vector<std::string> GetArgs() {return fArgs;}

        /**
         * Get the output of the last call to this module.
         * 
         * @returns The output of the last call.
        */
        const std::vector<float>& GetOutput() const {return fOutput;}

        /**
         * Get the output shape of the last call to this module.
         * 
         * @returns The output shape of the last call.
        */
        std::vector<int> GetShape() {return fOutShape;}

        virtual void Forward() = 0;  // Forward method to be implemented by each module.

        virtual std::vector<int> InferShape() = 0;  // Output shape inference to be implemented by each module.

        virtual std::string_view GetOperation() = 0;  // Operation name getter to be implemented by each module.
        
        virtual void SaveParameters(std::string dir) = 0;  // Parameter saver to be implemented by each module.

        virtual void LoadParameters() = 0;  // Parameter loader to be implemented by each module.

        virtual void LoadParameters(std::map<std::string, std::vector<float>>) = 0;  // Parameter loader to be implemented by each module.
    protected:
        std::string fName;  // Module name.  
        std::vector<std::shared_ptr<RModule>> fInputModules;  // Vector of input modules.
        std::vector<std::string> fInputs;  // Input names.
        std::vector<std::string> fArgs;  // Other arguments.
        std::vector<float> fOutput;  // Output of last call.
        std::vector<int> fOutShape;  // Output shape of last call.   
};

}  // TMVA.
}  // Experimental.
}  // SOFIE.

#endif  // TMVA_SOFIE_RMODULE_H_
