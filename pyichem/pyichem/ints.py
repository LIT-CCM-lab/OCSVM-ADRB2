import pandas as pd
import numpy as np
from grakel import Graph as grakelGraph
from pyichem.base_models import BatchCalculation
from biopandas.mol2 import PandasMol2
from scipy.spatial.distance import squareform, pdist
import os.path
import decimal
from itertools import permutations, combinations

import pdb

INTERACTIONS_PATH='ichem_outputs/interactions'

subst = {'CENT': ['SEC1', 'ALC2', 'ASC4', 'GLC5', 'PHC6'],
        'LIG': ['SEL1', 'ALL2', 'ASL4', 'GLL5', 'PHL6'],
        'PROT': ['SEP1', 'ALP2', 'ASP4', 'GLP5', 'PHP6'],
        'ELEC': ['SEC1', 'ALC2', 'LYC3', 'ASC4', 'PHC6', 'SEL1', 'ALL2', 'LYL3', 'ASL4', 'PHL6', 'SEP1', 'ALP2','LYP3', 'ASP4', 'PHP6'],
        'Pharma': ['SEL1', 'ALL2', 'LYL3', 'ASL4', 'PHL6']}

hb_simplification = {'SEC1': 'HBC8', 'ALC2': 'HBC8',
                    'SEP1': 'HBP8', 'ALP2': 'HBP8',
                    'SEL1': 'HBL8', 'ALL2': 'HBL8'}


class Ints(BatchCalculation):
    '''
    Class calling the Ints tool of IChem.

    WARNING !!! Metal chelation is not implemented in graph generation

    The supported definiton of IPAs during detection are available:
    MERG, IPAs placed on ligand, protein, and center
    CENT, IPAs placed in the middle point between the ligand interacting atom and the protein interacting atom
    PROT, IPAs placed on the protein interacting atom
    LIG, IPAs placed on the ligand iteracting atom

    From the IPAs different different interaction subgraphs can be generating, excluding some of the detected IPAs:
    CENT, IPAs placed in the middle point between the ligand interacting atom and the protein interacting atom
    PROT, IPAs placed on the protein interacting atom
    LIG, IPAs placed on the ligand iteracting atom
    ELEC, IPAs representing only polar interactions (hydrogen bonds, ionic bonds, pi-pi interactions)
    Pharma, IPAs placed on the ligand iteracting atom representing only polar interactions (hydrogen bonds, ionic bonds, pi-pi interactions)

    :param receptor_mol2: files containing the structure of the protein to use in IChem calculations
    :type receptor_mol2: list of str, optional
    :param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
    :type ligand_mol2: list of str, optional
    :param type_int: Type of protein-ligand interaction to detect.
    :type type_int: str, optional
    :param new_hyd: Use the NewHyd definition of hydrophobic contacts rather than the default one.
    :type new_hyd: bool, optional
    '''

    def __init__(self, receptor_mol2, ligand_mol2, type_int = 'MERG', new_hyd = False):
        if type_int in {'MERG', 'CENT', 'LIG', 'PROT'}:
            self.type_int = type_int
        else:
            raise ValueError(f'Interaction tpye {type_int} not recognized')
        if new_hyd:
            opt = f'-type {type_int} --newH'
        else:
            opt = f'-type {type_int}'

        super().__init__(receptor_mol2, ligand_mol2, INTERACTIONS_PATH, f'{INTERACTIONS_PATH}/out_ints_', 'ints', opt = opt)

    def compute_graphs(self, graph_type = 'grakel', threshold = None, subgraph = None, round_val = 1, simplify = False):
        '''
        Computes interaction graphs.
        Currently the output is given only as grakel graphs.

        :param graph_type: graph format used for the outputs, currently only grakel graphs are generated
        :type graph_type: str, optional
        :param threshold: distance threshold for edge deinition, currently not implemented in this version
        :type threshold: float, optional
        :param subgraph: extract a specific subgraph from the IPA data
        :type subgraph: bool
        :param round_val: closest value to which the Euclidean distance is approximated to define the edge weight
        :type round_val: float
        :param simplify: simplify the description of hydrogen bonds by treating hydrogen bond donor and acceptors as the same interaction type.
        :type simplify: bool
        :returns: Array containing the generated graphs
        :rtype: numpy array
        '''
        graphs = list()
        if protein:
            if len(self.output_location) != len(self.receptor_mol2):
                raise ValueError('The same number of interaction files and receptor structure file are required for residue relabelling')
            for file, prot_file in zip(self.output_location, self.receptor_mol2):
                graphs.append(graph_generator(file+f'_INTS_{self.type_int[0]}.mol2', graph_type = graph_type, threshold = threshold, subgraph = subgraph, round_val = round_val, simplify = simplify, protein_file = prot_file, backbone = backbone, mode = mode))
        else:
            for file in self.output_location:
                graphs.append(graph_generator(file+f'_INTS_{self.type_int[0]}.mol2', graph_type = graph_type, threshold = threshold, subgraph = subgraph, round_val = round_val, simplify = simplify, mode = mode))

        return np.array(graphs)

def graph_generator(file, mode = 'node', threshold = None, subgraph = None, round_val = None, simplify = False, graph_type = 'grakel'):
    '''
    Generates an interaction graph from a single .mol2 containing IPAs.
    Three different type of graphs are currently supported:
    node, only atom labels are set and edges are weighted
    edge, only edge labels are set
    node_edge, both node and edge labels are set 

    :param file: .mol2 file containing the IPAs location and types
    :type file: str
    :param mode: Describe the position of the labels in the graph
    :type mode: str
    :param threshold: distance threshold for edge deinition, currently not implemented in this version
    :type threshold: float, optional 
    :param subgraph: extract a specific subgraph from the IPA data
    :type subgraph: bool
    :param round_val: closest value to which the Euclidean distance is approximated to define the edge weight
    :type round_val: float
    :param simplify: simplify the description of hydrogen bonds by treating hydrogen bond donor and acceptors as the same interaction type.
    :type simplify: bool
    :param graph_type: graph format used for the outputs, currently only grakel graphs are generated
    :type graph_type: str, optional
    :returns: An interaction graph with the desired characteristics
    :rtype: grakel Graph
    '''
    func = {'node': graph_node_labels, 'edge': graph_edge_labels, 'node_edge': graph_node_edge_labels}
    if mode in func:
        dist, labels = graph_reader(file, threshold, subgraph, round_val, simplify, protein_file, backbone)
        return func[mode](dist, labels, graph_type = 'grakel')
    else:
        raise ValueError(f'Graph label mode {mode} is not recognized')

def graph_reader(file, threshold = None, subgraph = None, round_val = None, simplify = False):
    '''
    Function converting the IPAs in the .mol2 file to a list of labels and a distance matrix.
    :param threshold: distance threshold for edge deinition, currently not implemented in this version
    :type threshold: float, optional 
    :param subgraph: extract a specific subgraph from the IPA data
    :type subgraph: bool
    :param round_val: closest value to which the Euclidean distance is approximated to define the edge weight
    :type round_val: float
    :param simplify: simplify the description of hydrogen bonds by treating hydrogen bond donor and acceptors as the same interaction type.
    :type simplify: bool
    :returns: the upper triangular distance matrix as a 1D array, the list of node labels
    :rtype: numpy array, list
    '''
    if os.path.isfile(file):
        pmol = PandasMol2().read_mol2(file)
        if subgraph in subst:
            pmol._df = pmol.df[pmol.df['subst_name'].isin(subst[subgraph])]
        elif subgraph is not None:
            raise ValueError(f'The subgraph type {subgraph} is not supported')

        if simplify and pmol.df.shape[0] > 0:
            pmol._df = pmol.df.replace(hb_simplification)


        dist = pdist(pmol.df[['x', 'y', 'z']])
        labels_list = pmol.df['subst_name'].tolist()


        if round_val is not None:
            dist = np.rint(dist/round_val)*round_val
            dist[dist == 0] = round_val*0.1
     
    else:
        raise Exception(f'File {file} does not exist, impossible to compute the graph.')

    return dist, labels_list

def graph_node_labels(dist, node_labels, graph_type = 'grakel'):
    '''
    Function generating the interaction graph with labels on the nodes
    :param dist: triangular upper distance matrix as a 1D array
    :type dist: array
    :param node_labels: labels of the nodes in the graph
    :type node_labels: list
    :param graph_type: graph format used for the outputs, currently only grakel graphs are generated
    :type graph_type: str, optional
    :returns: Interaction graph
    :rtype: grakel Graph
    '''

    if graph_type == 'grakel':
    
        labels = {i:lab for i,lab in enumerate(node_labels)}
        sq_dist = squareform(dist)

        return grakelGraph(sq_dist, node_labels=labels)
    else:
        raise ValueError(f'Unsupported graph type {graph_type}')

def graph_edge_labels(dist, node_labels, graph_type = 'grakel'):
    '''
    Function generating the interaction graph with labels on the nodes.
    Each edge is described by a single label based on the two connected nodes and their distance.
    :param dist: triangular upper distance matrix as a 1D array
    :type dist: array
    :param node_labels: labels of the nodes in the graph
    :type node_labels: list
    :param graph_type: graph format used for the outputs, currently only grakel graphs are generated
    :type graph_type: str, optional
    :returns: Interaction graph
    :rtype: grakel Graph
    '''

    priority = {'SEC1':0, 'ALC2':1, 'LYC3':2, 'ASC4':3, 'GLC5':4, 'PHC6':5, 'SEL1':6, 'ALL2':7, 'LYL3':8, 'ASL4':9, 'GLL5':10, 'PHL6': 11, 'SEP1':12, 'ALP2':13,'LYP3':14, 'ASP4':15, 'GLP5':16, 'PHP6':17}

    if graph_type == 'grakel':
        
        permutation = combinations(node_labels, 2)
        idxs = combinations(np.arange(len(node_labels)), 2)
        adj = squareform(np.ones(len(dist)))

        edge_labels = dict()

        for node_labels, d, idx in zip(permutation, dist, idxs):
            edge_labels[(idx[0], idx[1])] = node_labels[0]+' '+node_labels[1]+' '+str(d)
            edge_labels[(idx[1], idx[0])] = node_labels[1]+' '+node_labels[0]+' '+str(d)

        return grakelGraph(adj, edge_labels = edge_labels)
    else:
        raise ValueError(f'Unsupported graph type {graph_type}')

def graph_node_edge_labels(dist, node_labels, graph_type = 'grakel'):
    '''
    Function generating the interaction graph with labels on the nodes
    An edge connecting node A and B is described by two labels:
    A B distance
    B A distance
    :param dist: triangular upper distance matrix as a 1D array
    :type dist: array
    :param node_labels: labels of the nodes in the graph
    :type node_labels: list
    :param graph_type: graph format used for the outputs, currently only grakel graphs are generated
    :type graph_type: str, optional
    :returns: Interaction graph
    :rtype: grakel Graph
    '''
    if graph_type == 'grakel':
    
        n_labels = {i:lab for i,lab in enumerate(node_labels)}
        sq_dist = squareform(dist)

        e_labels = {ind: sq_dist[ind[0], ind[1]] for ind in combinations(range(len(n_labels)), r = 2)}

        return grakelGraph(sq_dist, node_labels = n_labels, edge_labels = e_labels)
    else:
        raise ValueError(f'Unsupported graph type {graph_type}')