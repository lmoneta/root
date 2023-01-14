import ROOT 


def ParseFromMemory(GraphModule, GraphData):
    gin = ROOT.TMVA.Experimental.SOFIE.GNN_Init()
    gin.num_nodes = len(GraphData['nodes'])

    # extracting the edges
    edges = []
    for i in range(len(GraphData['senders'])):
        edges.append([GraphData['senders'][i],GraphData['receivers'][i]])
    gin.edges = edges

    gin.num_node_features = len(GraphData['nodes'][0])
    gin.num_edge_features = len(GraphData['edges'][0])
    gin.num_global_features = len(GraphData['globals'])

    gin.filename = GraphModule.name 

    # adding the node update function
    node_model = GraphModule._node_block._node_model
    if (node_model.name == "mlp"):
        num_layers = len(node_model._layers)
        upd = ROOT.TMVA.Experimental.SOFIE.RFunction_MLP(ROOT.TMVA.Experimental.SOFIE.FunctionTarget.NODES, num_layers, 0)
        kernel_tensor_names = []
        bias_tensor_names   = []

        for i in range(0, len(num_layers), 2):
            bias_tensor_names.append(node_model.variables[i])
            kernel_tensor_names.append(node_model.variables[i+1])
        upd.AddInitializedTensors([kernel_tensor_names, bias_tensor_names])
        gin.createUpdateFunction(upd)
    else:
        print("Invalid Model for node update.")
        return
    
    # adding the edge update function
    edge_model = GraphModule._edge_block._edge_model
    if (edge_model.name == "mlp"):
        num_layers = len(edge_model_layers)
        upd = ROOT.TMVA.Experimental.SOFIE.RFunction_MLP(ROOT.TMVA.Experimental.SOFIE.FunctionTarget.EDGES, num_layers, 0)
        kernel_tensor_names = []
        bias_tensor_names   = []

        for i in range(0, len(num_layers), 2):
            bias_tensor_names.append(edge_model.variables[i])
            kernel_tensor_names.append(edge_model.variables[i+1])
        upd.AddInitializedTensors([kernel_tensor_names, bias_tensor_names])
        gin.createUpdateFunction(upd)
    else:
        print("Invalid Model for edge update.")
        return

    # adding the global update function
    global_model = GraphModule._global_block._global_model
    if (global_model.name == "mlp"):
        num_layers = len(edge_model_layers)
        upd = ROOT.TMVA.Experimental.SOFIE.RFunction_MLP(ROOT.TMVA.Experimental.SOFIE.FunctionTarget.GLOBALS, num_layers, 0)
        kernel_tensor_names = []
        bias_tensor_names   = []

        for i in range(0, len(num_layers), 2):
            bias_tensor_names.append(global_model.variables[i])
            kernel_tensor_names.append(global_model.variables[i+1])
        upd.AddInitializedTensors([kernel_tensor_names, bias_tensor_names])
        gin.createUpdateFunction(upd)
    else:
        print("Invalid Model for edge update.")
        return
    
    # adding edge-node aggregate function
    edge_node_reducer = GraphModule._node_block._received_edges_aggregator._reducer.__qualname__
    if(edge_node_reducer == "unsorted_segment_sum"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_SUM()
    elif(edge_node_reducer == "unsorted_segment_mean"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_MEAN()
    else:
        print("Invalid aggregate function for edge-node reduction")
        return
    gin.createAggregateFunction(agg, ROOT.TMVA.Experimental.SOFIE.FunctionRelation.NODES_EDGES)

    
    # adding node-global aggregate function
    node_global_reducer = GraphModule._global_block._nodes_aggregator._reducer.__qualname__
    if(node_global_reducer == "unsorted_segment_sum"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_SUM()
    elif(node_global_reducer == "unsorted_segment_mean"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_MEAN()
    else:
        print("Invalid aggregate function for node-global reduction")
        return
    gin.createAggregateFunction(agg, ROOT.TMVA.Experimental.SOFIE.FunctionRelation.NODES_GLOBALS)

    # adding edge-global aggregate function
    node_global_reducer = GraphModule._global_block._edges_aggregator._reducer.__qualname__
    if(node_global_reducer == "unsorted_segment_sum"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_SUM()
    elif(node_global_reducer == "unsorted_segment_mean"):
        agg = ROOT.TMVA.Experimental.SOFIE.RFunction_MEAN()
    else:
        print("Invalid aggregate function for node-global reduction")
        return
    gin.createAggregateFunction(agg, ROOT.TMVA.Experimental.SOFIE.FunctionRelation.EDGES_GLOBALS)


    gnn_model = ROOT.TMVA.Experimental.SOFIE.RModel_GNN(gin)
    gnn_model.AddAddBlasRoutines(["gemm"])
    return gnn_model

    

