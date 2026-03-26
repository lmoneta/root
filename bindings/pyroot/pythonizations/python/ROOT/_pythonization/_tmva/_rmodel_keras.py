from tensorflow import keras
import os
import ROOT
import numpy as np
import math
import time
from .. import pythonization
from cppyy import gbl as gbl_namespace

def MakeKerasIdentity(layer):
    input = layer['layerInput']
    output = layer['layerOutput']
    fLayerType = layer_data['layerDType']
    fLayerInputName = input[0]
    fLayerOutputName = output[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Identity('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Identity does not yet support input type " + fLayerDType
        )

def MakeKerasBinary(layer):
    input = layer['layerInput']
    output = layer['layerOutput']
    fLayerType = layer_data['layerType'] 
    fLayerDType = layer_data['layerDType'] 
    fX1 = input[0]
    fX2 = input[1]
    fY = output[0]
    op = None
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        if fLayerType == "Add":
          op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_BasicBinary('Add')(fX1, fX2, fY)
        elif fLayerType == "Subtract":
          op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_BasicBinary('Sub')(fX1, fX2, fY)
        else:
          op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_BasicBinary('Mul')(fX1, fX2, fY)
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Identity does not yet support input type " + fLayerDType
        )
    return op


def MakeKerasConcat(layer):
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    attributes = layer['layerAttributes']
    input = [str(i) for i in finput]
    output = str(foutput[0])
    axis = int(attributes["axis"])
    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Concat('float')(inputs, axis, 0,  output)
    return op

def MakeKerasReshape(layer):
    """
    Create a Keras-compatible reshaping operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible reshaping operation using the SOFIE framework. Assumes layerDtype is float.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  name, data type, and other relevant information.

    Returns:
    ROperator_Reshape: A SOFIE framework operator representing the reshaping operation.
    """
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    attributes = layer['layerAttributes']
    flayername = attributes['_name']
    fOpMode =gbl_namespace.TMVA.Experimental.SOFIE.ReshapeOpMode.Reshape
    fLayerDType = layer['layerDType']
    fNameData = finput[0]
    fNameOutput = foutput[0]
    fNameShape = flayername + "ReshapeAxes"
    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Reshape('float')(fOpMode, 0, fNameData, fNameShape, fNameOutput)
    return op

def MakeKerasFlatten(layer):
    """
    Create a Keras-compatible flattening operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible flattening operation using the SOFIE framework.
    Flattening is the process of converting a multi-dimensional tensor into a
    one-dimensional tensor. Assumes layerDtype is float.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                name, data type, and other relevant information.

    Returns:
    ROperator_Reshape: A SOFIE framework operator representing the flattening operation.
    """
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    attributes = layer['layerAttributes']
    flayername = attributes['_name']
    fOpMode =gbl_namespace.TMVA.Experimental.SOFIE.ReshapeOpMode.Flatten
    fLayerDType = layer['layerDType']
    fNameData = finput[0]
    fNameOutput = foutput[0]
    fNameShape = flayername + "ReshapeAxes"
    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Reshape('float')(fOpMode, 0, fNameData, fNameShape, fNameOutput)
    return op


def MakeKerasBatchNorm(layer): 
    """
    Create a Keras-compatible batch normalization operation using SOFIE framework.

    This function takes a dictionary representing a batch normalization layer and its
    attributes and constructs a Keras-compatible batch normalization operation using
    the SOFIE framework. Batch normalization is used to normalize the activations of
    a neural network, typically applied after the convolutional or dense layers.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  gamma, beta, moving mean, moving variance, epsilon,
                  momentum, data type (assumed to be float), and other relevant information.

    Returns:
    ROperator_BatchNormalization: A SOFIE framework operator representing the batch normalization operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    attributes = layer['layerAttributes']
    gamma = attributes["gamma"]
    beta = attributes["beta"]
    moving_mean = attributes["moving_mean"]
    moving_variance = attributes["moving_variance"]
    fLayerDType = layer["layerDType"]
    fNX = str(finput[0])
    fNY = str(foutput[0])
    fNScale = str(gamma.name)
    fNB = str(beta.name)
    fNMean = str(moving_mean.name)
    fNVar = str(moving_variance.name)
    epsilon = attributes["epsilon"]
    momentum = attributes["momentum"]
    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_BatchNormalization('float')(epsilon, momentum, 0, fNX, fNScale, fNB, fNMean, fNVar, fNY)
    return op

def MakeKerasActivation(layer):
    attributes = layer['layerAttributes']
    activation = attributes['activation']
    fLayerActivation = str(activation.__name__)
    if fLayerActivation in mapKerasLayer.keys():
        return mapKerasLayer[fLayerActivation](layer)
    else:
        raise Exception("TMVA.SOFIE - parsing keras activation layer " + fLayerActivation + " is not yet supported")

def MakeKerasReLU(layer):
    """
    Create a Keras-compatible rectified linear unit (ReLU) activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible ReLU activation operation using the SOFIE framework.
    ReLU is a popular activation function that replaces all negative values in a tensor
    with zero, while leaving positive values unchanged.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type, which must be float.

    Returns:
    ROperator_Relu: A SOFIE framework operator representing the ReLU activation operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Relu('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Relu does not yet support input type " + fLayerDType
        )


def MakeKerasSeLU(layer):
    """
    Create a Keras-compatible scaled exponential linear unit (SeLU) activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible SeLU activation operation using the SOFIE framework.
    SeLU is a type of activation function that introduces self-normalizing properties
    to the neural network.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type - must be float32.

    Returns:
    ROperator_Selu: A SOFIE framework operator representing the SeLU activation operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Selu('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Selu does not yet support input type " + fLayerDType
        )


def MakeKerasSigmoid(layer):
    """
    Create a Keras-compatible sigmoid activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible sigmoid activation operation using the SOFIE framework.
    Sigmoid is a commonly used activation function that maps input values to the range
    between 0 and 1, providing a way to introduce non-linearity in neural networks.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type - must be float.

    Returns:
    ROperator_Sigmoid: A SOFIE framework operator representing the sigmoid activation operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Sigmoid('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Sigmoid does not yet support input type " + fLayerDType
        )


def MakeKerasSoftmax(layer):
    """
    Create a Keras-compatible softmax activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible softmax activation operation using the SOFIE framework.
    Softmax is an activation function that converts input values into a probability
    distribution, often used in the output layer of a neural network for multi-class
    classification tasks.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type - must be float.

    Returns:
    ROperator_Softmax: A SOFIE framework operator representing the softmax activation operation.
    """
    
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Softmax('float')(-1, fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Softmax does not yet support input type " + fLayerDType
        )


def MakeKerasLeakyRelu(layer):
    """
    Create a Keras-compatible Leaky ReLU activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible Leaky ReLU activation operation using the SOFIE framework.
    Leaky ReLU is a variation of the ReLU activation function that allows small negative
    values to pass through, introducing non-linearity while preventing "dying" neurons.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  attributes, and data type - must be float.

    Returns:
    ROperator_LeakyRelu: A SOFIE framework operator representing the Leaky ReLU activation operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    attributes = layer['layerAttributes']
    fAlpha = float(attributes["alpha"])
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_LeakyRelu('float')(fAlpha, fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator LeakyRelu does not yet support input type " + fLayerDType
        )


def MakeKerasTanh(layer):
    """
    Create a Keras-compatible hyperbolic tangent (tanh) activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible tanh activation operation using the SOFIE framework.
    Tanh is an activation function that squashes input values to the range between -1 and 1,
    introducing non-linearity in neural networks.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type - must be float.

    Returns:
    ROperator_Tanh: A SOFIE framework operator representing the tanh activation operation.
    """
        
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Tanh('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Tanh does not yet support input type " + fLayerDType
        )


def MakeKerasSwish(layer):
    """
    Create a Keras-compatible swish activation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible swish activation operation using the SOFIE framework.
    Swish is an activation function that aims to combine the benefits of ReLU and sigmoid,
    allowing some non-linearity while still keeping positive values unbounded.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  and data type.

    Returns:
    ROperator_Swish: A SOFIE framework operator representing the swish activation operation.
    """
    
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Swish('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Swish does not yet support input type " + fLayerDType
        )


def MakeKerasPermute(layer):
    """
    Create a Keras-compatible permutation operation using SOFIE framework.

    This function takes a dictionary representing a layer and its attributes and
    constructs a Keras-compatible permutation operation using the SOFIE framework.
    Permutation is an operation that rearranges the dimensions of a tensor based on
    specified dimensions.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  attributes, and data type - must be float.

    Returns:
    ROperator_Transpose: A SOFIE framework operator representing the permutation operation.
    """
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    attributes = layer['layerAttributes']
    fAttributePermute = np.asarray(attributes["dims"])
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        if len(fAttributePermute) > 0:
            op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')(fPermuteDims, fLayerInputName, fLayerOutputName)
        else:    
            op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')(fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Transpose does not yet support input type " + fLayerDType
        )


def MakeKerasDense(layer):
    """
    Create a Keras-compatible dense (fully connected) layer operation using SOFIE framework.

    This function takes a dictionary representing a dense layer and its attributes and
    constructs a Keras-compatible dense (fully connected) layer operation using the SOFIE framework.
    A dense layer applies a matrix multiplication between the input tensor and weight matrix,
    and adds a bias term.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  layer weight names, and data type - must be float.

    Returns:
    ROperator_Gemm: A SOFIE framework operator representing the dense layer operation.
    """  
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    fWeightNames = layer["layerWeight"]
    fKernelName = fWeightNames[0]
    fBiasName = fWeightNames[1]
    attr_alpha = 1.0
    attr_beta  = 1.0
    attr_transA = 0
    attr_transB = 0
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Gemm['float'](attr_alpha, attr_beta, attr_transA, attr_transB, fLayerInputName, fKernelName, fBiasName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Gemm does not yet support input type " + fLayerDType
        )


def MakeKerasConv(layer): 
    """
    Create a Keras-compatible convolutional layer operation using SOFIE framework.

    This function takes a dictionary representing a convolutional layer and its attributes and
    constructs a Keras-compatible convolutional layer operation using the SOFIE framework.
    A convolutional layer applies a convolution operation between the input tensor and a set
    of learnable filters (kernels).

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  data type (must be float), weight and bias name, kernel size, dilations, padding and strides. 
                  When padding is same (keep in the same dimensions), the padding shape is calculated.

    Returns:
    ROperator_Conv: A SOFIE framework operator representing the convolutional layer operation.
    """
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerDType = layer['layerDType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    attributes = layer['layerAttributes']
    fWeightNames = layer["layerWeight"]
    fKernelName = fWeightNames[0]
    fBiasName = fWeightNames[1]
    fAttrDilations = attributes["dilation_rate"]
    fAttrGroup = int(attributes["groups"])
    fAttrKernelShape = attributes["kernel_size"]
    fKerasPadding = str(attributes["padding"])
    fAttrStrides = attributes["strides"]
    
    if fKerasPadding == 'valid':
        fAttrAutopad = 'VALID'
    elif fKerasPadding == 'same':
        fAttrAutopad = 'NOTSET'
        fInputShape = attributes['_build_input_shape']
        inputHeight = fInputShape[1]
        inputWidth = fInputShape[2]
        outputHeight = math.ceil(float(inputHeight) / float(fAttrStrides[0]))
        outputWidth = math.ceil(float(inputWidth) / float(fAttrStrides[1]))
        padding_height = max((outputHeight - 1) * fAttrStrides[0] + fAttrKernelShape[0] - inputHeight, 0)
        padding_width = max((outputWidth - 1) * fAttrStrides[1] + fAttrKernelShape[1] - inputWidth, 0)
        padding_top = math.floor(padding_height / 2)
        padding_bottom = padding_height - padding_top
        padding_left = math.floor(padding_width / 2)
        padding_right = padding_width - padding_left
        fAttrPads = [padding_top, padding_bottom, padding_left, padding_right]
    else:
        raise RuntimeError(
            "TMVA::SOFIE - RModel Keras Parser doesn't yet supports Convolution layer with padding " + fKerasPadding
        )
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Conv['float'](fAttrAutopad, fAttrDilations, fAttrGroup, 
                                                                  fAttrKernelShape, fAttrPads, fAttrStrides, 
                                                                  fLayerInputName, fKernelName, fBiasName, 
                                                                  fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Gemm does not yet support input type " + fLayerDType
        )


def MakeKerasPooling(layer):
    """
    Create a Keras-compatible pooling layer operation using SOFIE framework.

    This function takes a dictionary representing a pooling layer and its attributes and
    constructs a Keras-compatible pooling layer operation using the SOFIE framework.
    Pooling layers downsample the input tensor by selecting a representative value from
    a group of neighboring values, either by taking the maximum or the average.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  layer type (the selection rule), the pool size, padding, strides, and data type.

    Returns:
    ROperator_Pool: A SOFIE framework operator representing the pooling layer operation.
    """
    
    #extract attributes from layer data
    fLayerDType = layer['layerDType']
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    fLayerType = layer['layerType']
    fLayerInputName = finput[0]
    fLayerOutputName = foutput[0]
    pool_atrr =gbl_namespace.TMVA.Experimental.SOFIE.RAttributes_Pool()
    attributes = layer['layerAttributes']
    fAttrKernelShape = attributes["pool_size"]
    fKerasPadding = str(attributes["padding"])
    fAttrStrides = attributes["strides"]
    if fKerasPadding == 'valid':
        fAttrAutopad = 'VALID'
    elif fKerasPadding == 'same':
        fAttrAutopad = 'NOTSET'
    else:
        raise RuntimeError(
            "TMVA::SOFIE - RModel Keras Parser doesn't yet supports Convolution layer with padding " + fKerasPadding
        )
    pool_atrr.dilations = list(fAttrDilations)
    pool_atrr.strides = list(fAttrStrides)
    pool_atrr.pads = fpads
    pool_atrr.kernel_shape = list(fAttrKernelShape)
    pool_atrr.auto_pad = fAttrAutopad    
    
    #choose pooling type
    if fLayerType.startswith("Max"):
        PoolMode = gbl_namespace.TMVA.Experimental.SOFIE.PoolOpMode.MaxPool
    elif fLayerType.startswith("AveragePool"):
        PoolMode = gbl_namespace.TMVA.Experimental.SOFIE.PoolOpMode.AveragePool
    elif fLayerType.startswith("GlobalAverage"):
        PoolMode = gbl_namespace.TMVA.Experimental.SOFIE.PoolOpMode.GloabalAveragePool
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator poolong does not yet support pooling type " + fLayerType
        )
    
    #Set default values
    fAttrDilations = (1,1)
    fpads = [0,0,0,0,0,0]
    pool_atrr.ceil_mode = 0
    pool_atrr.count_include_pad = 0
    pool_atrr.storage_order = 0
    
    #create operator
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Pool['float'](PoolMode, pool_atrr, fLayerInputName, fLayerOutputName)
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator Pooling does not yet support input type " + fLayerDType
        )


def MakeKerasRNN(layer): 
    """
    Create a Keras-compatible RNN (Recurrent Neural Network) layer operation using SOFIE framework.

    This function takes a dictionary representing an RNN layer and its attributes and
    constructs a Keras-compatible RNN layer operation using the SOFIE framework.
    RNN layers are used to model sequences, and they maintain internal states that are
    updated through recurrent connections.

    Parameters:
    layer (dict): A dictionary containing layer information including input, output,
                  layer type, attributes, weights, and data type - must be float.

    Returns:
    ROperator_RNN: A SOFIE framework operator representing the RNN layer operation.
    """
    
    # Extract required information from the layer dictionary
    fLayerDType = layer['layerDType']
    finput = layer['layerInput']
    foutput = layer['layerOutput']
    attributes = layer['layerAttributes']
    direction = attributes['direction']
    hidden_size = attributes["hidden_size"]
    layout = int(attributes["layout"])
    nameX = finput[0]
    nameY = foutput[0]
    nameW = layer["layerWeight"][0]
    nameR = layer["layerWeight"][1]
    if len(layer["layerWeight"]) > 2:
        nameB = layer["layerWeight"][2]
    else:
        nameB = ""
    
    # Check if the provided activation function is supported
    fPActivation = attributes['activation']
    if not fPActivation.__name__ in ['relu', 'sigmoid', 'tanh', 'softsign', 'softplus']: #avoiding functions with parameters
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator RNN does not yet support activation function " + fPActivation.__name__
        )
    activations = [fPActivation.__name__[0].upper()+fPActivation.__name__[1:]]

    #set default values
    activation_alpha = {}
    activation_beta = {}
    clip = 0.0
    nameY_h = ""
    nameInitial_h = ""
    name_seq_len = ""
    
    if gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fLayerDType) == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
        if layer['layerType'] == "SimpleRNN":
            op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_RNN['float'](activation_alpha, activation_beta, activations, clip, direction, hidden_size, layout, nameX, nameW, nameR, nameB, name_seq_len, nameInitial_h, nameY, nameY_h)
        
        elif layer['layerType'] == "GRU":
            #an additional activation function is required, given by the user
            activations.insert(0,attributes['recurrent_activation'])
            
            #new variable needed:
            linear_before_reset = attributes['linear_before_reset']
            op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_GRU['float'](activation_alpha, activation_beta, activations, clip, direction, hidden_size, layout, linear_before_reset, nameX, nameW, nameR, nameB, name_seq_len, nameInitial_h, nameY, nameY_h)
        
        elif layer['layerType'] == "LSTM":
            #an additional activation function is required, the first given by the user, the second set to tanh as default
            fPRecurrentActivation = attributes['recurrent_activation']
            if not fPActivation.__name__ in ['relu', 'sigmoid', 'tanh', 'softsign', 'softplus']: #avoiding functions with parameters
                raise RuntimeError(
                    "TMVA::SOFIE - Unsupported - Operator RNN does not yet support recurrent activation function " + fPActivation.__name__
                )
            fPRecurrentActivationName = fPRecurrentActivation.__name__[0].upper()+fPRecurrentActivation.__name__[1:]
            activations.insert(0,fPRecurrentActivationName)
            activations.insert(2,'Tanh')            
            
            #new variables needed:
            input_forget = 0
            nameInitial_c = ""
            nameP = "" #No peephole connections in keras LSTM model
            nameY_c = ""
            op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_LSTM['float'](activation_alpha, activation_beta, activations, clip, direction, hidden_size, input_forget, layout, nameX, nameW, nameR, nameB, name_seq_len, nameInitial_h, nameInitial_c, nameP, nameY, nameY_h, nameY_c)
        
        else: 
            raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator RNN does not yet support operator type " + layer['layerType']
        ) 
        return op
    else:
        raise RuntimeError(
            "TMVA::SOFIE - Unsupported - Operator RNN does not yet support input type " + fLayerDType
        )   

#Set global dictionaries, mapping layers to corresponding functions that create their ROperator instances
mapKerasLayer = {"Activation": MakeKerasActivation,
                 "Permute": MakeKerasPermute,
                 "BatchNormalization": MakeKerasBatchNorm,
                 "Reshape": MakeKerasReshape,
                 "Flatten": MakeKerasFlatten,
                 "Concatenate": MakeKerasConcat,
                 "swish": MakeKerasSwish,
                 "Add": MakeKerasBinary,
                 "Subtract": MakeKerasBinary,
                 "Multiply": MakeKerasBinary,
                 "Softmax": MakeKerasSoftmax,
                 "tanh": MakeKerasTanh,
                 "Identity": MakeKerasIdentity,
                 "Dropout": MakeKerasIdentity,
                 "ReLU": MakeKerasReLU,
                 "relu": MakeKerasReLU,
                 "selu": MakeKerasSeLU,
                 "sigmoid": MakeKerasSigmoid,
                 "LeakyReLU": MakeKerasLeakyRelu, 
                 "softmax": MakeKerasSoftmax, 
                 "MaxPooling2D": MakeKerasPooling,
                 "SimpleRNN": MakeKerasRNN,
                 "GRU": MakeKerasRNN,
                 "LSTM": MakeKerasRNN,
                 }

mapKerasLayerWithActivation = {"Dense": MakeKerasDense,"Conv2D": MakeKerasConv}


def add_layer_into_RModel(rmodel, layer_data):
    """
    Add a Keras layer operation to an existing RModel using the SOFIE framework.

    This function takes an existing RModel and a dictionary representing a Keras layer
    and its attributes, and adds the corresponding layer operation to the RModel using
    the SOFIE framework. The function supports various types of Keras layers, including
    those with or without activation functions.

    Parameters:
    rmodel (RModel): An existing RModel to which the layer operation will be added.
    layer_data (dict): A dictionary containing layer information including type,
                      attributes, input, output, and layer data type.

    Returns:
    RModel: The updated RModel after adding the layer operation.

    Raises exception: If the provided layer type or activation function is not supported.
    """
    
    fLayerType = layer_data['layerType']
    
    #reshape and flatten layers don't have weights, but they are needed inside the list of initialized tensor list in the Rmodel
    if fLayerType == "Reshape" or fLayerType == "Flatten":
        Attributes = layer_data['layerAttributes']
        LayerName = Attributes['_name']
        if fLayerType == "Reshape":
            TargetShape = np.asarray(Attributes['target_shape']).astype("int")
            TargetShape = np.insert(TargetShape,0,0)
        else:
            input_shape = layer_data['layerAttributes']['_build_input_shape']
            TargetShape = [gbl_namespace.TMVA.Experimental.SOFIE.ConvertShapeToLength(input_shape[1:])]
            TargetShape = np.asarray(TargetShape)
        
        #since the AddInitializedTensor method in RModel requires unique pointer, we call a helper function in c++ that does the conversion from a regular pointer to unique one in c++
        rmodel.AddInitializedTensorFromPy['long'](LayerName+"ReshapeAxes",gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.INT64,[len(TargetShape)], TargetShape)
    
    #These layers only have one operator - excluding the recurrent layers, in which the activation function(s) are included in the recurrent operator
    if fLayerType in mapKerasLayer.keys():
        Attribues = layer_data['layerAttributes']
        inputs = layer_data['layerInput']
        outputs = layer_data['layerOutput']
        LayerName = Attribues['_name']
        
        #Pooling layers in keras by default assume the channels dimension is the last one, 
        #while in onnx (and the RModel) it is the first one (other than batch size), 
        #so a transpose is needed before and after the pooling, if the data format is channels last (can be set to channels first by the user).
        if fLayerType == 'MaxPooling2D':
            if layer_data['channels_last']:
                op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')([0,3,1,2], inputs[0], LayerName+"PreTrans")
                rmodel.AddOperatorFromPy(op)
                inputs[0] = LayerName+"PreTrans"
                layer_data["layerInput"] = inputs
                outputs[0] = LayerName+fLayerType
                layer_data['layerOutput'] = outputs
        rmodel.AddOperatorFromPy(mapKerasLayer[fLayerType](layer_data))
        if fLayerType == 'MaxPooling2D':
            if layer_data['channels_last']:
                op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')([0,2,3,1], LayerName+fLayerType, LayerName+"PostTrans")
                rmodel.AddOperatorFromPy(op)
        return rmodel
    
    #These layers require two operators - dense/conv and their activation funciton
    elif fLayerType in mapKerasLayerWithActivation.keys():
        Attribues = layer_data['layerAttributes']
        LayerName = Attribues['_name']
        fPActivation = Attribues['activation']
        LayerActivation = fPActivation.__name__
        if LayerActivation in ['selu', 'sigmoid']:
            rmodel.AddNeededStdLib("cmath")
        
        #if there is an activation function after the layer
        if LayerActivation != 'linear':
            outputs = layer_data['layerOutput']
            inputs = layer_data['layerInput']
            fActivationLayerOutput = outputs[0]
            
            #like pooling, convolutional layer from keras requires transpose before and after to match the onnx format 
            # if the data format is channels last (can be set to channels first by the user).
            if fLayerType == 'Conv2D':
                if layer_data['channels_last']:
                    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')([0,3,1,2], inputs[0], LayerName+"PreTrans")
                    rmodel.AddOperatorFromPy(op)
                    inputs[0] = LayerName+"PreTrans"
                    layer_data["layerInput"] = inputs
            outputs[0] = LayerName+fLayerType
            layer_data['layerOutput'] = outputs
            op = mapKerasLayerWithActivation[fLayerType](layer_data)
            rmodel.AddOperatorFromPy(op)
            Activation_layer_input = LayerName+fLayerType
            if fLayerType == 'Conv2D':
                if layer_data['channels_last']:
                    op = gbl_namespace.TMVA.Experimental.SOFIE.ROperator_Transpose('float')([0,2,3,1], LayerName+fLayerType, LayerName+"PostTrans")
                    rmodel.AddOperatorFromPy(op)
                    Activation_layer_input = LayerName + "PostTrans"
            
            #Adding the activation function
            inputs[0] = Activation_layer_input
            outputs[0] = fActivationLayerOutput
            layer_data['layerInput'] = inputs
            layer_data['layerOutput'] = outputs
            if not LayerActivation in mapKerasLayer.keys():
                raise Exception("TMVA.SOFIE - parsing keras activation function " + LayerActivation + " is not yet supported")
            rmodel.AddOperatorFromPy(mapKerasLayer[LayerActivation](layer_data))
            
        else: #there is a bug here if it is conv and the activation is linear, need to add transpose before and after
            rmodel.AddOperatorFromPy(mapKerasLayerWithActivation[fLayerType](layer_data))
        return rmodel
    else:
        raise Exception("TMVA.SOFIE - parsing keras layer " + fLayerType + " is not yet supported")


def Keras_Parser_into_RModel(filename):
    #Check if file exists
    if not os.path.exists(filename):
        raise RuntimeError("Model file {} not found!".format(filename))
        
    #load model
    keras_model = keras.models.load_model(modelFile)
    keras_model.load_weights(modelFile)
    
    #create new RModel object
    sep = '/'
    if os.name == 'nt':
        sep = '\\'
    
    isep = filename.rfind(sep)
    filename_nodir = filename
    if isep != -1:
        filename_nodir = filename[isep+1:]
    
    ttime = time.time()
    gmt_time = time.gmtime(ttime)
    parsetime = time.asctime(gmt_time)
    
    rmodel = gbl_namespace.TMVA.Experimental.SOFIE.RModel.RModel(filename_nodir, parsetime)
    
    #iterate over the layers and add them to the RModel
    for layer in keras_model.layers:
        layer_data={}
        layer_data['layerType']=layer.__class__.__name__
        layer_data['layerAttributes']=layer.__dict__
        layer_data['layerInput']=[x.name for x in layer.input] if isinstance(layer.input,list) else [layer.input.name]
        layer_data['layerOutput']=[x.name for x in layer.output] if isinstance(layer.output,list) else [layer.output.name]
        layer_data['layerDType']=layer.dtype
        layer_data['layerWeight']=[x.name for x in layer.weights]
        
        #for convolutional and pooling layers we need to know the format of the data
        if layer_data['layerType'] in ['Conv2D', 'MaxPooling2D']:
            layer_data['channels_last'] = True if layer.data_format == 'channels_last' else False
            
        #for recurrent type layers we need to extract additional unique information
        if layer_data['layerType'] in ["SimpleRNN", "LSTM", "GRU"]:
            layer_data['layerAttributes']['activation'] = layer.activation
            layer_data['layerAttributes']['direction'] = 'backward' if layer.go_backwards else 'forward'
            layer_data['layerAttributes']["units"] = layer.units
            layer_data['layerAttributes']["layout"] = layer.input.shape[0] is None
            layer_data['layerAttributes']["hidden_size"] = layer.output.shape[-1]
            
            #for GRU and LSTM we need to extract an additional activation function
            if layer_data['layerType'] != "SimpleRNN": 
                layer_data['layerAttributes']['recurrent_activation'] = layer.recurrent_activation
            
            #for GRU there are two variants of the reset gate location, we need to know which one is it
            if layer_data['layerType'] == "GRU":
                layer_data['layerAttributes']['linear_before_reset'] = 1 if layer.reset_after and layer.recurrent_activation.__name__ == "sigmoid" else 0
                       
        if layer_data['layerInput'][0].startswith('max_pooling2d'):
            pooling_layer_name = layer_data['layerInput'][0].split('/')[0]
            layer_data['layerInput'][0] = pooling_layer_name + 'PostTrans'
        
        fLayerType = layer_data['layerType']
        #Ignoring the input layer for models built using Keras Functional API
        #NEED TO TEST KERAS FUNCTIONAL API
        if(fLayerType == "InputLayer"):
            continue;

        #Adding any required routines depending on the Layer types for generating inference code.
        elif (fLayerType == "Dense"):
            rmodel.AddBlasRoutines({"Gemm", "Gemv"})
        elif (fLayerType == "BatchNormalization"):
            rmodel.AddBlasRoutines({"Copy", "Axpy"})
        elif (fLayerType == "Conv1D" or fLayerType == "Conv2D" or fLayerType == "Conv3D"):
            rmodel.AddBlasRoutines({"Gemm", "Axpy"})
        rmodel = add_layer_into_RModel(rmodel, layer_data)

    # Extracting model's weights
    weight = []
    for idx in range(len(keras_model.get_weights())):
        weightProp = {}
        weightProp['name'] = keras_model.weights[idx].name
        weightProp['dtype'] = keras_model.get_weights()[idx].dtype.name
        if 'conv' in keras_model.weights[idx].name and keras_model.weights[idx].shape.ndims == 4:
            weightProp['value'] = keras_model.get_weights()[idx].transpose((3, 2, 0, 1)).copy()
        else:
            weightProp['value'] = keras_model.get_weights()[idx]
        weight.append(weightProp)

    # Traversing through all the Weight tensors
    for weightIter in range(len(weight)):
        fWeightTensor = weight[weightIter]
        fWeightName = fWeightTensor['name']
        fWeightDType =gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fWeightTensor['dtype'])
        fWeightTensorValue = fWeightTensor['value']
        fWeightTensorSize = 1
        fWeightTensorShape = []
        
        #IS IT BATCH SIZE? CHECK ONNX
        if fWeightName.startswith("simple_rnn") or fWeightName.startswith("lstm") or (fWeightName.startswith("gru") and not 'bias' in fWeightName):
            fWeightTensorShape.append(1)
        
        # Building the shape vector and finding the tensor size
        for j in range(len(fWeightTensorValue.shape)):
            fWeightTensorShape.append(fWeightTensorValue.shape[j])
            fWeightTensorSize *= fWeightTensorValue.shape[j]
        
        if fWeightDType == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
            fWeightArray = fWeightTensorValue
            
            #weights conversion format between keras and onnx for lstm: the order of the different elements (input, output, forget, cell) inside the vector/matrix is different
            if fWeightName.startswith("lstm"):
                if 'kernel' in fWeightName:
                    units = int(fWeightArray.shape[1]/4)
                    W_i = fWeightArray[:, :units].copy()
                    W_f = fWeightArray[:, units: units * 2].copy()
                    W_c = fWeightArray[:, units * 2: units * 3].copy()
                    W_o = fWeightArray[:, units * 3:].copy()
                    fWeightArray[:, units: units * 2] = W_o
                    fWeightArray[:, units * 2: units * 3] = W_f
                    fWeightArray[:, units * 3:] = W_c
                else: #bias
                    units = int(fWeightArray.shape[0]/4)
                    W_i = fWeightArray[:units].copy()
                    W_f = fWeightArray[units: units * 2].copy()
                    W_c = fWeightArray[units * 2: units * 3].copy()
                    W_o = fWeightArray[units * 3:].copy()
                    fWeightArray[units: units * 2] = W_o
                    fWeightArray[units * 2: units * 3] = W_f
                    fWeightArray[units * 3:] = W_c
           
            #need to make specific adjustments for recurrent weights and biases
            if (fWeightName.startswith("simple_rnn") or fWeightName.startswith("lstm") or fWeightName.startswith("gru")):
                #reshaping weight matrices for recurrent layers due to keras-onnx inconsistencies
                if 'kernel' in fWeightName:
                    fWeightArray = np.transpose(fWeightArray)
                    fWeightTensorShape[1], fWeightTensorShape[2] = fWeightTensorShape[2], fWeightTensorShape[1]
                
                fData = fWeightArray.flatten()
                
                #the recurrent bias and the cell bias can be the same, in which case we need to add a vector of zeros for the recurrent bias
                if 'bias' in fWeightName and len(fData.shape) == 1:
                    fWeightTensorShape[1] *= 2
                    fRbias = fData.copy()*0
                    fData = np.concatenate((fData,fRbias))

            else:
                fData = fWeightArray.flatten()
                
            rmodel.AddInitializedTensorFromPy['float'](fWeightName, gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT, fWeightTensorShape, fData)
        else:
            raise TypeError("Type error: TMVA SOFIE does not yet support data layer type: " + fWeightDType)
    
    # Extracting input tensor info
    fPInputs = keras_model.input_names
    fPInputShape = keras_model.input_shape if isinstance(keras_model.input_shape, list) else [keras_model.input_shape]
    fPInputDType = []
    for idx in range(len(keras_model.inputs)):
        fPInputDType.append(keras_model.inputs[idx].dtype.__str__()[9:-2])
    
    if len(fPInputShape) == 1:
        fInputName = fPInputs[0]
        fInputDType =gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fPInputDType[0])
        if fInputDType == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
            if fPInputShape[0][0] is None or fPInputShape[0][0] <= 0:
                fPInputShape = list(fPInputShape[0])
                fPInputShape[0] = 1
            rmodel.AddInputTensorInfo(fInputName, gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT, fPInputShape)
            rmodel.AddInputTensorName(fInputName) 
        else:
            raise TypeError("Type error: TMVA SOFIE does not yet support data type "+TMVA.Experimental.SOFIE.ConvertStringToType(fInputDType))
    else:
         #Iterating through multiple input tensors
        for fInputName, fInputDType, fInputShapeTuple in zip(fPInputs, fPInputDType, fPInputShape):
            fInputDType =gbl_namespace.TMVA.Experimental.SOFIE.ConvertStringToType(fInputDType)
            if fInputDType == gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT:
                if fInputShapeTuple[0] is None or fInputShapeTuple[0] <= 0:
                    fInputShapeTuple = list(fInputShapeTuple)
                    fInputShapeTuple[0] = 1
                    print("Model does not have a defined batch size. Assuming it is 1 - input shape: ", fInputShapeTuple)
                rmodel.AddInputTensorInfo(fInputName, gbl_namespace.TMVA.Experimental.SOFIE.ETensorType.FLOAT, fInputShapeTuple)
                rmodel.AddInputTensorName(fInputName)
            else:
                raise TypeError("Type error: TMVA SOFIE does not yet support data type "+TMVA.Experimental.SOFIE.ConvertStringToType(fInputDType))             
        
    # Adding OutputTensorInfos
    outputNames = []
    for layerName in keras_model.output_names:
        outputNames.append(keras_model.get_layer(layerName).output.name)
    rmodel.AddOutputTensorNameList(outputNames)
    return rmodel

@pythonization("Keras_Parser_into_RModel", ns="TMVA::Experimental")
def Keras_Parser_into_RModel(klass):
    # Parameters:
    # klass: class to be pythonized
    klass.__init__ = Keras_Parser_into_RModel

