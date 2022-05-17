import pandas as pd
import numpy as np
from grakel import Graph as grakelGraph
from pyichem.base_models import BatchCalculation
from biopandas.mol2 import PandasMol2
from scipy.spatial.distance import squareform, pdist
import os.path

INTERACTIONS_PATH='ichem_outputs/interactions'

class Ints(BatchCalculation):
	'''
	Class calling the Ints tool of IChem

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
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
		graphs = list()
		for file in self.output_location:
			graphs.append(self.__graph(file+f'_INTS_{self.type_int[0]}.mol2', graph_type = graph_type, threshold = threshold, subgraph = subgraph, round_val = round_val, simplify = simplify))

		return np.array(graphs)

	def __graph(self, file, graph_type = 'grakel', threshold = None, subgraph = None, round_val = 1, simplify = False):
		subst = {'CENT': ['SEC1', 'ALC2', 'ASC4', 'GLC5', 'PHC6'],
		'LIG': ['SEL1', 'ALL2', 'ASL4', 'GLL5', 'PHL6'],
		'PROT': ['SEP1', 'ALP2', 'ASP4', 'GLP5', 'PHP6'],
		'PC': ['SEP1', 'ALP2', 'ASP4', 'GLP5', 'PHP6', 'SEC1', 'ALC2', 'ASC4', 'GLC5', 'PHC6'],
		'ELEC': ['SEC1', 'ALC2', 'LYC3', 'ASC4', 'PHC6', 'SEL1', 'ALL2', 'LYL3', 'ASL4', 'PHL6', 'SEP1', 'ALP2','LYP3', 'ASP4', 'PHP6'],
		'EPL': ['SEL1', 'ALL2', 'ASL4', 'PHL6', 'SEP1', 'ALP2', 'ASP4', 'PHP6'],
		'ELIG': ['SEL1', 'ALL2', 'ASL4', 'PHL6'],
		'HYD': ['GLC5', 'GLP5', 'GLL5'],
		'HB': ['SEC1', 'ALC2', 'SEL1', 'ALL2', 'SEP1', 'ALP2'],
		'CHARGE': ['ASC4', 'ASL4', 'ASP4','LYC3', 'LYL3', 'LYP3'],
		'ARO': ['PHC6', 'PHL6', 'PHP6']}
		if os.path.isfile(file):
			pmol = PandasMol2().read_mol2(file)
			if subgraph in subst:
				pmol._df = pmol.df[pmol.df['subst_name'].isin(subst[subgraph])]
			elif subgraph is not None:
				raise ValueError(f'The subgraph type {subgraph} is not supported')
			dist = pdist(pmol.df[['x', 'y', 'z']])
			dist = np.rint(dist/round_val)*round_val
			dist[dist == 0] = round_val*0.1
			dist = squareform(dist)

			
			if graph_type == 'grakel':
				if simplify:
					labels_list = np.array(pmol.df['subst_name'].tolist())
					if len(labels_list) > 0:
						labels_list[(labels_list == 'SEC1') | (labels_list == 'ALC2')] = 'HBC8'
						labels_list[(labels_list == 'SEP1') | (labels_list == 'ALP2')] = 'HBP8'
						labels_list[(labels_list == 'SEL1') | (labels_list == 'ALL2')] = 'HBL8'
						labels_list[(labels_list == 'ASC4') ] = 'INC9'
						labels_list[(labels_list == 'ASP4') ] = 'INP9'
						labels_list[(labels_list == 'ASL4') ] = 'INL9'
					labels = {i:lab for i,lab in enumerate(labels_list)}

				else:
					labels = {i:lab for i,lab in enumerate(pmol.df['subst_name'].tolist())}
				return grakelGraph(dist, node_labels=labels)
			else:
				raise ValueError(f'Unsupported graph type {graph_type}')
		else:
			raise Exception(f'File {file} does not exist, impossible to compute the graph.')