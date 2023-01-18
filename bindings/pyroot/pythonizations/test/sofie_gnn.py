import unittest
import ROOT

import numpy as np
import graph_nets as gn
from graph_nets import utils_tf
import sonnet as snt


GLOBAL_FEATURE_SIZE = 2
NODE_FEATURE_SIZE = 2
EDGE_FEATURE_SIZE = 2


# for generating input data for graph nets, 
# from https://github.com/deepmind/graph_nets/blob/master/graph_nets/demos/graph_nets_basics.ipynb
def get_graph_data_dict(num_nodes, num_edges):
  return {
      "globals": np.random.rand(GLOBAL_FEATURE_SIZE).astype(np.float32),
      "nodes": np.random.rand(num_nodes, NODE_FEATURE_SIZE).astype(np.float32),
      "edges": np.random.rand(num_edges, EDGE_FEATURE_SIZE).astype(np.float32),
      "senders": np.random.randint(num_nodes, size=num_edges, dtype=np.int32),
      "receivers": np.random.randint(num_nodes, size=num_edges, dtype=np.int32),
  }


class SOFIE_GNN(unittest.TestCase):
    """
    Tests for the pythonizations of ParseFromMemory method of SOFIE GNN. 
    """

    def test_parse_gnn(self):
        '''
        Test that parsed GNN model from a graphnets model generates correct 
        inference code
        '''
        GraphModule = gn.modules.GraphNetwork(
            edge_model_fn=lambda: snt.nets.MLP([2,2], activate_final=True),
            node_model_fn=lambda: snt.nets.MLP([2,2], activate_final=True),
            global_model_fn=lambda: snt.nets.MLP([2,2], activate_final=True))

        GraphData = get_graph_data_dict(2,1)
        input_graphs = utils_tf.data_dicts_to_graphs_tuple([GraphData])
        output = GraphModule(input_graphs)
        
        # Parsing model to RModel_GNN
        model = ROOT.TMVA.Experimental.SOFIE.RModel_GNN.ParseFromMemory(GraphModule, GraphData)
        model.Generate()
        model.OutputGenerated()

        ROOT.gInterpreter.Declare('#include "graph_network.hxx"')
        input_data = ROOT.TMVA.Experimental.SOFIE.GNN_Data()

        for i in GraphData['nodes'].flatten():
            input_data.node_data.push_back(i)

        for i in GraphData['edges'].flatten():
            input_data.edge_data.push_back(i)

        for i in GraphData['globals'].flatten():
            input_data.global_data.push_back(i)

        ROOT.TMVA_SOFIE_graph_network.infer(input_data)
        
        output_node_data = output.nodes.numpy().flatten()
        output_edge_data = output.edges.numpy().flatten()
        output_global_data = output.globals.numpy().flatten()

        for i in range(len(output_node_data)):
            self.assertAlmostEqual(output_node_data[i], input_data.node_data[i], 2)

        for i in range(len(output_edge_data)):
            self.assertAlmostEqual(output_edge_data[i], input_data.edge_data[i], 2)

        for i in range(len(output_global_data)):
            self.assertAlmostEqual(output_global_data[i], input_data.global_data[i], 2)




if __name__ == '__main__':
    unittest.main()
