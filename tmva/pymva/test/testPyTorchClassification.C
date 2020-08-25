#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"

TString pythonSrc = "import torch\n\
from torch import nn\n\
\n\
# Define model\n\
model = nn.Sequential(\n\
                nn.Linear(4, 64),\n\
                nn.ReLU(),\n\
                nn.Linear(64, 2),\n\
                nn.Softmax(dim=1))\n\
\n\
# Construct loss function and Optimizer.\n\
criterion = torch.nn.MSELoss()\n\
optimizer = torch.optim.SGD\n\
\n\
\n\
def fit(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):\n\
    trainer = optimizer(model.parameters(), lr=0.01)\n\
    schedule, schedulerSteps = scheduler\n\
    best_val = None\n\
\n\
    for epoch in range(num_epochs):\n\
        # Training Loop\n\
        # Set to train mode\n\
        model.train()\n\
        running_train_loss = 0.0\n\
        running_val_loss = 0.0\n\
        for i, (X, y) in enumerate(train_loader):\n\
            trainer.zero_grad()\n\
            output = model(X)\n\
            train_loss = criterion(output, y)\n\
            train_loss.backward()\n\
            trainer.step()\n\
\n\
            # print train statistics\n\
            running_train_loss += train_loss.item()\n\
            if i % 32 == 31:    # print every 32 mini-batches\n\
                print(f\"[{epoch+1}, {i+1}] train loss: {running_train_loss / 32 :.3f}\")\n\
                running_train_loss = 0.0\n\
\n\
        if schedule:\n\
            schedule(optimizer, epoch, schedulerSteps)\n\
\n\
        # Validation Loop\n\
        # Set to eval mode\n\
        model.eval()\n\
        with torch.no_grad():\n\
            for i, (X, y) in enumerate(val_loader):\n\
                output = model(X)\n\
                val_loss = criterion(output, y)\n\
                running_val_loss += val_loss.item()\n\
\n\
            curr_val = running_val_loss / len(val_loader)\n\
            if save_best:\n\
               if best_val==None:\n\
                   best_val = curr_val\n\
               best_val = save_best(model, curr_val, best_val)\n\
\n\
            # print val statistics per epoch\n\
            print(f\"[{epoch+1}] val loss: {curr_val :.3f}\")\n\
            running_val_loss = 0.0\n\
\n\
    print(f\"Finished Training on {epoch+1} Epochs!\")\n\
\n\
    return model\n\
\n\
\n\
def predict(model, test_X, batch_size=32):\n\
    # Set to eval mode\n\
    model.eval()\n\
   \n\
    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))\n\
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)\n\
\n\
    predictions = []\n\
    with torch.no_grad():\n\
        for i, data in enumerate(test_loader):\n\
            X = data[0]\n\
            outputs = model(X)\n\
            predictions.append(outputs)\n\
        preds = torch.cat(predictions)\n\
   \n\
    return preds.numpy()\n\
\n\
\n\
load_model_custom_objects = {\"optimizer\": optimizer, \"criterion\": criterion, \"train_func\": fit, \"predict_func\": predict}\n\
\n\
# Store model to file\n\
m = torch.jit.script(model)\n\
torch.jit.save(m,\"PyTorchModelClassification.pt\")\n";


int testPyTorchClassification(){
   // Get data file
   std::cout << "Get test data..." << std::endl;
   TString fname = "./tmva_class_example.root";
   if (gSystem->AccessPathName(fname))  // file does not exist in local directory
      gSystem->Exec("curl -O http://root.cern.ch/files/tmva_class_example.root");
   TFile *input = TFile::Open(fname);

   // Build model from python file
   std::cout << "Generate PyTorch model..." << std::endl;
   UInt_t ret;
   ret = gSystem->Exec("echo '"+pythonSrc+"' > generatePyTorchModelClassification.py");
   if(ret!=0){
       std::cout << "[ERROR] Failed to write python code to file" << std::endl;
       return 1;
   }
   ret = gSystem->Exec("python generatePyTorchModelClassification.py");
   if(ret!=0){
       std::cout << "[ERROR] Failed to generate model using python" << std::endl;
       return 1;
   }

   // // Setup PyMVA and factory
   std::cout << "Setup TMVA..." << std::endl;
   TMVA::PyMethodBase::PyInitialize();
   TFile* outputFile = TFile::Open("ResultsTestPyTorchClassification.root", "RECREATE");
   TMVA::Factory *factory = new TMVA::Factory("testPyTorchClassification", outputFile,
      "!V:Silent:Color:!DrawProgressBar:AnalysisType=Classification");

   // Load data
   TMVA::DataLoader *dataloader = new TMVA::DataLoader("datasetTestPyTorchClassification");

   TTree *signal = (TTree*)input->Get("TreeS");
   TTree *background = (TTree*)input->Get("TreeB");
   dataloader->AddSignalTree(signal);
   dataloader->AddBackgroundTree(background);

   dataloader->AddVariable("var1");
   dataloader->AddVariable("var2");
   dataloader->AddVariable("var3");
   dataloader->AddVariable("var4");

   dataloader->PrepareTrainingAndTestTree("",
      "SplitMode=Random:NormMode=NumEvents:!V");

   // Book and train method
   factory->BookMethod(dataloader, TMVA::Types::kPyTorch, "PyTorch",
      "!H:!V:VarTransform=D,G:FilenameModel=PyTorchModelClassification.pt:FilenameTrainedModel=trainedPyTorchModelClassification.pt:NumEpochs=10:BatchSize=32:UserCode=generatePyTorchModelClassification.py");
   std::cout << "Training model..." << std::endl;
   factory->TrainAllMethods();

   // Clean-up
   delete factory;
   delete dataloader;
   delete outputFile;

   // Setup reader
   UInt_t numEvents = 100;
   std::cout << "Run reader and classify " << numEvents << " events..." << std::endl;
   TMVA::Reader *reader = new TMVA::Reader("!Color:Silent");
   Float_t vars[4];
   reader->AddVariable("var1", vars+0);
   reader->AddVariable("var2", vars+1);
   reader->AddVariable("var3", vars+2);
   reader->AddVariable("var4", vars+3);
   reader->BookMVA("PyTorch", "datasetTestPyTorchClassification/weights/testPyTorchClassification_PyTorch.weights.xml");

   // Get mean response of method on signal and background events
   signal->SetBranchAddress("var1", vars+0);
   signal->SetBranchAddress("var2", vars+1);
   signal->SetBranchAddress("var3", vars+2);
   signal->SetBranchAddress("var4", vars+3);

   background->SetBranchAddress("var1", vars+0);
   background->SetBranchAddress("var2", vars+1);
   background->SetBranchAddress("var3", vars+2);
   background->SetBranchAddress("var4", vars+3);

   Float_t meanMvaSignal = 0;
   Float_t meanMvaBackground = 0;
   for(UInt_t i=0; i<numEvents; i++){
      signal->GetEntry(i);
      meanMvaSignal += reader->EvaluateMVA("PyTorch");
      background->GetEntry(i);
      meanMvaBackground += reader->EvaluateMVA("PyTorch");
   }
   meanMvaSignal = meanMvaSignal/float(numEvents);
   meanMvaBackground = meanMvaBackground/float(numEvents);

   // Check whether the response is obviously better than guessing
   std::cout << "Mean MVA response on signal: " << meanMvaSignal << std::endl;
   if(meanMvaSignal < 0.6){
      std::cout << "[ERROR] Mean response on signal is " << meanMvaSignal << " (<0.6)" << std::endl;
      return 1;
   }
   std::cout << "Mean MVA response on background: " << meanMvaBackground << std::endl;
   if(meanMvaBackground > 0.4){
      std::cout << "[ERROR] Mean response on background is " << meanMvaBackground << " (>0.4)" << std::endl;
      return 1;
   }

   return 0;
}

int main(){
   int err = testPyTorchClassification();
   return err;
}
