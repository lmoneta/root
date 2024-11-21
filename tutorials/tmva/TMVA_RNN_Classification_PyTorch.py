import ROOT
from ROOT import TMVA

ninput = 30

ntime = 10

batchSize = 100
maxepochs = 20

use_type = 1

nTotEvts = 10000 # total events to be generated for signal or background

useKeras = True


useTMVA_RNN = True
useTMVA_DNN = True
useTMVA_BDT = False

rnn_types = ["RNN", "LSTM", "GRU"]
use_rnn_type = [1, 1, 1]
if (use_type >=0 & use_type < 3) :
      use_rnn_type = [0,0,0]
      use_rnn_type[use_type] = 1;

archString = "CPU"
writeOutputFile = True

rnn_type = "RNN"

ROOT.TMVA.Tools.Instance()
ROOT.TMVA.PyMethodBase.PyInitialize()

num_threads = 0   # use by default all threads
#    do enable MT running
if (num_threads >= 0):
    ROOT.EnableImplicitMT(num_threads)
    if (num_threads > 0):
        ROOT.gSystem.Setenv("OMP_NUM_THREADS", num_threads)
    else:
      ROOT.gSystem.Setenv("OMP_NUM_THREADS", "1")


print("Running with nthreads  = " + str(ROOT.GetThreadPoolSize()) + "\n" )

inputFileName = "time_data_t10_d30.root"

fileExist = ROOT.gSystem.AccessPathName(inputFileName)

#if file does not exists create it
if (fileExist==None):
    MakeTimeData(nTotEvts,ntime, ninput)

inputFile = ROOT.TFile.Open(inputFileName)
if (inputFile==None):
    Error("TMVA_RNN_Classification", "Error opening input file %s - exit", inputFileName.Data())

print("--- RNNClassification  : Using input file: " + inputFile.GetName()+"\n")

#   Create a ROOT output file where TMVA will store ntuples, histograms, etc.
outfileName = "data_RNN_"+ archString +".root"

if (writeOutputFile):
    outputFile = ROOT.TFile.Open(outfileName, "RECREATE")

#  Creating the factory object
factory = ROOT.TMVA.Factory("TMVAClassification", outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=None:!Correlations:"+"AnalysisType=Classification:ModelPersistence")

dataloader =ROOT.TMVA.DataLoader("dataset")

signalTree = inputFile.Get("sgn")
background = inputFile.Get("bkg")

signalTree.Print()
nvar = ninput * ntime

# add variables - use new AddVariablesArray function
for i in range(ntime):
    varName = "vars_time"+str(i)
    dataloader.AddVariablesArray(varName,ninput,'F')

dataloader.AddSignalTree(signalTree, 1.0)
dataloader.AddBackgroundTree(background, 1.0)

# check given input
datainfo = dataloader.GetDataSetInfo()
vars = datainfo.GetListOfVariables()
print("number of variables is " + str(vars.size())+ "\n")
for v in vars:
    print(str(v)+"\n")

nTrainSig = 0.8 * nTotEvts
nTrainBkg = 0.8 *  nTotEvts

#build the string options for DataLoader::PrepareTrainingAndTestTree
prepareOptions = "nTrain_Signal="+str(nTrainSig)+":nTrain_Background="+str(nTrainBkg)+":SplitMode=Random:SplitSeed=100:NormMode=NumEvents:!V:!CalcCorrelations"


# Apply additional cuts on the signal and background samples (can be different)
mycuts = ROOT.TCut("")   ## for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
mycutb = ROOT.TCut("")   ## for example: TCut mycutb = "abs(var1)<0.5";

dataloader.PrepareTrainingAndTestTree(mycuts, mycutb, prepareOptions)

print("prepared DATA LOADER " )

if (useTMVA_RNN):
    for i in range(3):
        if (use_rnn_type[i]==None):
            continue
        rnn_type = str(rnn_types[i])

#          define the inputlayout string for RNN
#          the input data should be organize as   following:
#          input layout for RNN:    time x ndim

        inputLayoutString = "InputLayout="+str(ntime)+"|"+str(ninput)

        # Define RNN layer layout
        # it should be   LayerType (RNN or LSTM or GRU) |  number of units | number of inputs | time steps | remember output (typically no=0 | return full sequence
        rnnLayout = str(rnn_type) + "|10|"+ str(ninput) + "|" + str(ntime) + "|0|1"

        #        add after RNN a reshape layer (needed top flatten the output) and a dense layer with 64 units and a last one
        #        Note the last layer is linear because  when using Crossentropy a Sigmoid is applied already
        layoutString ="Layout=" + rnnLayout + ",RESHAPE|FLAT,DENSE|64|TANH,LINEAR"

        #Defining Training strategies. Different training strings can be concatenate. Use however only one
        trainingString1 = "LearningRate=1e-3,Momentum=0.0,Repetitions=1,"+"ConvergenceSteps=5,BatchSize="+str(batchSize)+",TestRepetitions=1,"+"WeightDecay=1e-2,Regularization=None,MaxEpochs="+str(maxepochs
        )+","+"Optimizer=ADAM,DropConfig=0.0+0.+0.+0."

        trainingStrategyString="TrainingStrategy="
        trainingStrategyString += trainingString1; # + "|" + trainingString2

        # Define the full RNN Noption string adding the final options for all network
        rnnOptions = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"+"WeightInitialization=XAVIERUNIFORM:ValidationSize=0.2:RandomSeed=1234"
        rnnOptions +=  ":" + inputLayoutString
        rnnOptions +=  ":" + layoutString
        rnnOptions +=  ":" + trainingStrategyString
        rnnOptions +=  ":" + "Architecture=" + str(archString)

        rnnName = "TMVA_" + rnn_type
        factory.BookMethod(dataloader, TMVA.Types.kDL, rnnName, rnnOptions)

if (useTMVA_DNN):
#    Method DL with Dense Layer
    inputLayoutString = "InputLayout=1|1|" + str(ntime * ninput)

    layoutString = "Layout=DENSE|64|TANH,DENSE|TANH|64,DENSE|TANH|64,LINEAR"
#   Training strategies.
    trainingString1 = "LearningRate=1e-3,Momentum=0.0,Repetitions=1,"+"ConvergenceSteps=10,BatchSize=256,TestRepetitions=1,"+"WeightDecay=1e-4,Regularization=None,MaxEpochs=20"+"DropConfig=0.0+0.+0.+0.,Optimizer=ADAM"
    trainingStrategyString = "TrainingStrategy="
    trainingStrategyString += trainingString1 # + "|" + trainingString2

      # General Options.
    dnnOptions = "!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=None:"+"WeightInitialization=XAVIER:RandomSeed=0" 

    dnnOptions +=  ":" + inputLayoutString
    dnnOptions +=  ":" + layoutString
    dnnOptions +=  ":" + trainingStrategyString
    dnnOptions +=  ":" + "Architecture=" + str(archString)


    dnnName = "TMVA_DNN"
    factory.BookMethod(dataloader, TMVA.Types.kDL, dnnName, dnnOptions)
   

import torch
from torch import nn

# Define model

# Custom Reshape Layer
class Reshape(torch.nn.Module):
    def forward(self, x):
        return x.view(-1,1,10,30)

# CNN Model Definition
net = torch.nn.Sequential(
    Reshape(),
    nn.Conv2d(1, 10, kernel_size=3, padding=1),
    nn.ReLU(),
    nn.BatchNorm2d(10),
    nn.Conv2d(10, 10, kernel_size=3, padding=1),
    nn.ReLU(),
    nn.MaxPool2d(kernel_size=2),
    nn.Flatten(),
    nn.Linear(10*5*15, 300),
    nn.ReLU(),
    nn.Linear(300, 2),
    nn.Sigmoid()
    )

# Construct loss function and Optimizer.
criterion = nn.BCELoss()
optimizer = torch.optim.Adam


def fit(model, train_loader, val_loader, num_epochs, batch_size, optimizer, criterion, save_best, scheduler):
    trainer = optimizer(model.parameters(), lr=0.01)
    schedule, schedulerSteps = scheduler
    best_val = None

    # Setup GPU
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    for epoch in range(num_epochs):
        # Training Loop
        # Set to train mode
        model.train()
        running_train_loss = 0.0
        running_val_loss = 0.0
        for i, (X, y) in enumerate(train_loader):
            trainer.zero_grad()
            X, y = X.to(device), y.to(device)
            output = model(X)
            target = y
            train_loss = criterion(output, target)
            train_loss.backward()
            trainer.step()

            # print train statistics
            running_train_loss += train_loss.item()
            if i % 4 == 3:    # print every 4 mini-batches
                print(f"[{epoch+1}, {i+1}] train loss: {running_train_loss / 4 :.3f}")
                running_train_loss = 0.0

        if schedule:
            schedule(optimizer, epoch, schedulerSteps)

        # Validation Loop
        # Set to eval mode
        model.eval()
        with torch.no_grad():
            for i, (X, y) in enumerate(val_loader):
                X, y = X.to(device), y.to(device)
                output = model(X)
                target = y
                val_loss = criterion(output, target)
                running_val_loss += val_loss.item()

            curr_val = running_val_loss / len(val_loader)
            if save_best:
               if best_val==None:
                   best_val = curr_val
               best_val = save_best(model, curr_val, best_val)

            # print val statistics per epoch
            print(f"[{epoch+1}] val loss: {curr_val :.3f}")
            running_val_loss = 0.0

    print(f"Finished Training on {epoch+1} Epochs!")

    return model


def predict(model, test_X, batch_size=100):
    # Set to eval mode

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)

    model.eval()


    test_dataset = torch.utils.data.TensorDataset(torch.Tensor(test_X))
    test_loader = torch.utils.data.DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    predictions = []
    with torch.no_grad():
        for i, data in enumerate(test_loader):
            X = data[0].to(device)
            outputs = model(X)
            predictions.append(outputs)
        preds = torch.cat(predictions)

    return preds.cpu().numpy()


load_model_custom_objects = {"optimizer": optimizer, "criterion": criterion, "train_func": fit, "predict_func": predict}

# Store model to file
m = torch.jit.script(net)
torch.jit.save(m,"PyTorchModelRNN.pt")

factory.BookMethod(dataloader, ROOT.TMVA.Types.kPyTorch, "PyTorch","H:!V:VarTransform=None:FilenameModel=PyTorchModelRNN.pt:" + "FilenameTrainedModel=PyTorchTrainedModelRNN.pt:NumEpochs=20:BatchSize=100")    

# Train all methods
factory.TrainAllMethods()

print("nthreads  = "+ str(ROOT.GetThreadPoolSize()) + "\n")

# Evaluate all MVAs using the set of test events
factory.TestAllMethods()

# Evaluate and compare performance of all configured MVAs
factory.EvaluateAllMethods()
#  check method

#  plot ROC curve
c1 = factory.GetROCCurve(dataloader)
c1.Draw() 


if (outputFile):
    outputFile.Close()
