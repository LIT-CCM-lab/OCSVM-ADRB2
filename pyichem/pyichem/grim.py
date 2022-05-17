import pandas as pd
import numpy as np
from pyichem.base_models import BatchCalculation

GRIM_PATH = 'ichem_outputs/GRIM'

class GrimStructures(BatchCalculation):
	'''
	Class calling the Grim tool of IChem.
	

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
	:param ifp_format: format of the IFP
	:type ifp-format: str, optional
	:param output_file: name of the file where all the generated ifp are stored
	:type output_file: str
	'''

	def __init__(self, receptor_mol2, ligand_mol2, ref_receptor_mol2, ref_ligand_mol2, values = True, inter = 'MERG', match = 'MERG', score = 'FCT', ref_name = None, comp_name = None):
		'''Constructor method'''
		
		self.ref_receptor_mol2 = ref_receptor_mol2
		self.ref_ligand_mol2 = ref_ligand_mol2
		self.inter = inter
		self.match = match
		self.score = score
		self.values = values
		self.ref_name = ref_name
		self.comp_name = comp_name
		option = self.__check_option()
		super().__init__(receptor_mol2, ligand_mol2, GRIM_PATH, '', 'grim', output_f = False, opt = option)

	def _write_input(self, input_file):
		'''
		Write the input file given to IChem.

		:param input_file: name of the input file
		:type input_file: str
		'''

		with open(f'{self.folder}/{input_file}', 'w') as input_f:
			if self.ref_name is not None and self.comp_name is not None:
				for receptor, ligand, c_name in zip(self.receptor_mol2, self.ligand_mol2, self.comp_name):
					for ref_receptor, ref_ligand, r_name in zip(self.ref_receptor_mol2, self.ref_ligand_mol2, self.ref_name):
						r_name = r_name.replace('\n', '')
						c_name = c_name.replace('\n', '')
						input_f.write(f'{self.options} -rn {r_name} -cn {c_name} {self.software} {ref_receptor} {ref_ligand} {receptor} {ligand}\n')
			elif self.ref_name is not None:
				for receptor, ligand, c_name in zip(self.receptor_mol2, self.ligand_mol2, self.comp_name):
					for ref_receptor, ref_ligand in zip(self.ref_receptor_mol2, self.ref_ligand_mol2):
						c_name = c_name.replace('\n', '')
						input_f.write(f'{self.options} -rn {r_name} {self.software} {ref_receptor} {ref_ligand} {receptor} {ligand}\n')
			elif self.ref_name is not None:
				for receptor, ligand in zip(self.receptor_mol2, self.ligand_mol2):
					for ref_receptor, ref_ligand, r_name in zip(self.ref_receptor_mol2, self.ref_ligand_mol2, self.ref_name):
						r_name = r_name.replace('\n', '')
						input_f.write(f'{self.options} -cn {c_name} {self.software} {ref_receptor} {ref_ligand} {receptor} {ligand}\n')
			else:
				for receptor, ligand in zip(self.receptor_mol2, self.ligand_mol2):
					for ref_receptor, ref_ligand in zip(self.ref_receptor_mol2, self.ref_ligand_mol2):
						input_f.write(f'{self.options} {self.software} {ref_receptor} {ref_ligand} {receptor} {ligand}\n')


		self.input_file = input_file

	def map_results(self):
		'''
		Generates a pandas DataFrame containing and the protein and ligand files used to generate a specific row in the dataframe containing all the generated IFPs.

		:return: table containing the files used to generate an IFP and the index of the IFP in the final table
		:rtype: pandas DataFrame
		'''
		index_lines = np.arange(0, len(self.receptor_mol2)*len(self.ref_receptor_mol2))

		receptor_files = list()
		ligand_files = list()

		for _ in range(len(self.ref_receptor_mol2)):
			receptor_files = receptor_files + self.receptor_mol2
			ligand_files = ligand_files + self.ligand_mol2

		ref_receptor = list()
		ref_ligand = list()
		for ref_l, ref_r in zip(self.ref_ligand_mol2, self.ref_receptor_mol2):
			for _ in range(len(self.receptor_mol2)):
				ref_receptor.append(ref_r)
				ref_ligand.append(ref_l)

		try:
			return pd.DataFrame(data = {'Receptor_file' : receptor_files, 'Ligand_file' : ligand_files, 'Receptor_file_reference' : ref_receptor, 'Ligand_file_reference' : ref_ligand, 'Grim_index' : index_lines})
		except ValueError:
			raise ValueError("Arrays must all be same length \n Check IChem output file to detect missing results \n If there are missing values try to run the calculation manually using IChem with the generated input file \n Possible segmentation fault for one of the calculations")

	def __check_option(self):
		option = ''
		allowed_values = ['MERG', 'LIG', 'CENT', 'PROT']

		if self.values:
			option = option + '--values '
		elif self.inter != 'MERG':
			if self.inter in allowed_values:
				option = option + f'-outInt {self.inter} '
			else:
				raise ValueError(f'Unrecognized value for the output: {self.inter}')
		if self.match != 'MERG':
			if self.match in allowed_values:
				option = option + f'-match {self.match} '
			else:
				raise ValueError(f'Unrecognized value for the alignment: {self.align}')
		if self.score != 'FCT':
			if self.score == 'STD':
				option = option + '-score STD '
			else:
				raise ValueError(f'Unrecognized value for the score: {self.score}')

		return option


class GrimIPA(BatchCalculation):

	'''
	Class calling the Grim tool of IChem.
	

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
	:param ifp_format: format of the IFP
	:type ifp-format: str, optional
	:param output_file: name of the file where all the generated ifp are stored
	:type output_file: str
	'''

	def __init__(self, comp_ipa, ref_ipa, values = True, inter = 'MERG', match = 'MERG', score = 'FCT', ref_name = None, comp_name = None):
		'''Constructor method'''
		
		self.comp_ipa = comp_ipa
		self.ref_ipa = ref_ipa
		self.inter = inter
		self.match = match
		self.score = score
		self.values = values
		self.ref_name = ref_name
		self.comp_name = comp_name
		option = self.__check_option()
		super().__init__([], [], GRIM_PATH, '', 'grim', output_f = False, opt = option)

	def _write_input(self, input_file):
		'''
		Write the input file given to IChem.

		:param input_file: name of the input file
		:type input_file: str
		'''

		with open(f'{self.folder}/{input_file}', 'w') as input_f:
			if self.ref_name is not None and self.comp_name is not None:
				for ipa, c_name in zip(self.comp_ipa, self.comp_name):
					for ripa, r_name in zip(self.ref_ipa, self.ref_name):
						r_name = r_name.replace('\n', '')
						c_name = c_name.replace('\n', '')
						input_f.write(f'{self.options} -rn {r_name} -cn {c_name} {self.software} {ripa} {ipa}\n')
			elif self.comp_name is not None:
				for ipa, c_name in zip(self.comp_ipa, self.comp_name):
					for ripa in self.ref_ipa:
						c_name = c_name.replace('\n', '')
						input_f.write(f'{self.options} -cn {c_name} {self.software} {ripa} {ipa}\n')
			elif self.ref_name is not None:
				for ipa in self.comp_ipa:
					for ripa, r_name in zip(self.ref_ipa, self.ref_name):
						r_name = r_name.replace('\n', '')
						input_f.write(f'{self.options} -rn {r_name} {self.software} {ripa} {ipa}\n')
			else:
				for ipa in self.comp_ipa:
					for ripa in self.ref_ipa:
						input_f.write(f'{self.options} {self.software} {ripa} {ipa}\n')


		self.input_file = input_file

	def map_results(self):
		'''
		Generates a pandas DataFrame containing and the protein and ligand files used to generate a specific row in the dataframe containing all the generated IFPs.

		:return: table containing the files used to generate an IFP and the index of the IFP in the final table
		:rtype: pandas DataFrame
		'''
		index_lines = np.arange(0, len(self.comp_ipa)*len(self.ref_ipa))

		ipa_files = list()

		for _ in range(len(self.ref_ipa)):
			ipa_files = ipa_files + self.comp_ipa

		ref_ipa_files = list()
		for ref_ipa in self.ref_ipa:
			for _ in range(len(self.comp_ipa)):
				ref_ipa_files.append(ref_ipa)

		try:
			return pd.DataFrame(data = {'Compound_IPA' : ipa_files, 'Reference_IPA' : ref_ipa_files, 'Grim_index' : index_lines})
		except ValueError:
			raise ValueError("Arrays must all be same length \n Check IChem output file to detect missing results \n If there are missing values try to run the calculation manually using IChem with the generated input file \n Possible segmentation fault for one of the calculations")

	def __check_option(self):
		option = ''
		allowed_values = ['MERG', 'LIG', 'CENT', 'PROT']

		if self.values:
			option = option + '--values '
		elif self.inter != 'MERG':
			if self.inter in allowed_values:
				option = option + f'-outInt {self.inter} '
			else:
				raise ValueError(f'Unrecognized value for the output: {self.inter}')
		if self.match != 'MERG':
			if self.match in allowed_values:
				option = option + f'-match {self.match} '
			else:
				raise ValueError(f'Unrecognized value for the alignment: {self.align}')
		if self.score != 'FCT':
			if self.score == 'STD':
				option = option + '-score STD '
			else:
				raise ValueError(f'Unrecognized value for the score: {self.score}')

		return option