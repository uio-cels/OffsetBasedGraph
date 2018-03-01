from .graph import BaseGraph, Graph, BlockCollection, Block, BlockArray, AdjListAsNumpyArrays
import numpy as np
import json
from collections import defaultdict
import logging
import pickle

class GraphWithReversals(Graph):
    def __init__(self, blocks, adj_list,
                 create_reverse_adj_list=True,
                 rev_adj_list=None):
        super(GraphWithReversals, self).__init__(
            blocks, adj_list,
            create_reverse_adj_list=create_reverse_adj_list,
            rev_adj_list=rev_adj_list)
