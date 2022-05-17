import pandas as pd
import numpy as np
import subprocess
import sys
import pdb
from pyichem.base_models import BatchCalculation

IFP_PATH = 'ichem_outputs/IFP'

class Ifp(BatchCalculation):
	'''
	Class calling the Ifp tool of IChem.
	Different types of interactions can be included in the IFP, look at the IChem manual for more information.

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
	:param ifp_format: format of the IFP
	:type ifp-format: str, optional
	:param output_file: name of the file where all the generated ifp are stored
	:type output_file: str
	'''

	def __init__(self, receptor_mol2, ligand_mol2, ifp_format = 'regular', output_file = 'ligands.ifp'):
		'''Constructor method'''
		formats = {'regular': '', 'polar': '--polar', 'extended': '--extended'}
		self.ifp_format = ifp_format
		if ifp_format in formats:
			self.ifp_option = formats[ifp_format]
		else:
			raise ValueError(f'Invalid fingerprint type for {ifp_format}')
		super().__init__(receptor_mol2, ligand_mol2, IFP_PATH, f'{IFP_PATH}/ligands.ifp', 'IFP', output_f = False, stdout_name = f'{IFP_PATH}/{output_file}', opt = self.ifp_option)

	def read_ifp(self):
		'''
		Reads the generated IFP files in a single pandas DataFrame.
		The function reads both a series of file or a single file.
		If a file is missing or it is empty NaN values are given to all the interactions

		:returns: fingerprints attribute of the IFP object
		:rtype: pandas DataFrame
		'''

		if hasattr(self, 'output_location'):
			for file in self.output_location:
				if hasattr(self, 'fingerprints'):
					if file_check(file):
						self.fingerprints = self.fingerprints.append(ifp_reader(file, self.ifp_format), ignore_index = True)
					else:
						self.fingerprints = self.fingerprints.append(pd.DataFrame(np.full((1,self.fingerprints.shape[1]), np.nan), columns = self.fingerprints.columns), ignore_index = True)
				else:
					self.fingerprints = ifp_reader(file, self.ifp_format)
		else:
			self.fingerprints = ifp_reader(self.stdout, self.ifp_format)

		return self.fingerprints

	def map_results(self):
		'''
		Generates a pandas DataFrame containing and the protein and ligand files used to generate a specific row in the dataframe containing all the generated IFPs.

		:return: table containing the files used to generate an IFP and the index of the IFP in the final table
		:rtype: pandas DataFrame
		'''
		index_lines = np.arange(0, self.fingerprints.shape[0])
		try:
			return pd.DataFrame(data = {'Receptor_file' : self.receptor_mol2, 'Ligand_file' : self.ligand_mol2, 'IFP_index' : index_lines})
		except ValueError:
			raise ValueError("Arrays must all be same length \n Check IChem output file to detect missing results \n If there are missing values try to run the calculation manually using IChem with the generated input file \n Possible segmentation fault for one of the calculations")


	def fp_interaction(self, interactions):
		'''
		Filters the IFP table selecting only the specified interactions types.

		:param interactions: selected interaction(s)
		:type interactions: str or list of str
		:return: table containing the selected interactions
		:rtype: pandas DataFrame
		'''
		if not isinstance(interactions, list):
			interactions = [interactions]
		return filter_interaction(self.fingerprints, interactions)

	def fp_residues(self, residues):
		'''
		Filters the IFP table selecting only the specified residues.

		:param residues: selected residue(s)
		:type residues: str or list of str
		:return: table containing the selected residues
		:rtype: pandas DataFrame
		'''
		if not isinstance(residues, list):
			residues = [residues]
		return filter_residues(self.fingerprints, residues)

	def calculate_lbl(self, input_file = 'ichem_input.in', stdout_c = True, stderr_c = True):
		'''
		Generates the input file, then launch the calculations line by line, this is not the optimal
		way to perform calculations but it is a workaround for 'Segmentation fault (core dumped)'
		errors observed when launching from file.

		:param input_file: name of the input file
		:type input_file: str
		:param stdout_c: capture the stdout stream from IChem
		:type stdout_c: bool
		:param stderr_c: capture the stderr stream from IChem
		:type stderr_c: bool
		'''

		self._write_input(input_file)

		with open(f'{self.folder}/{self.input_file}', 'r') as inpt:
			for i,line in enumerate(inpt):
				command = line.split()
				command.insert(0, self.ichem_path)
				try:
					process_output = subprocess.run(command,
			 			stdout = subprocess.PIPE if stdout_c else None,
			 			stderr = subprocess.PIPE if stderr_c else None)
				except FileNotFoundError:
					raise FileNotFoundError(f'The selected IChem path was incorrect: {self.ichem_path}.\nPlease update the file at: {os.path.dirname(os.path.realpath(__file__))}/software_path.yml')


				
				if process_output.stdout is not None:
					with open(f'{IFP_PATH}/ligands.ifp', 'a') as out:
						out.write(f'\nIFP from {line}\n')
						msg = process_output.stdout.decode('utf-8')
						if msg != '':
							out.write(msg)
						else:
							out.write(f'|WARNING IChem was not able to calculate the IFP with the line : {line} \n')


				if process_output.stderr is not None:
					with open('ichem_stderr.txt', 'a') as err:
						err.write(process_output.stderr.decode('utf-8'))

def ifp_reader(file, ifp_type = 'regular'):
	'''
	Read an IFP file containing one or multiple IFPs

	:param file: file containing the IFPs
	:type file: str
	:param ifp_type: type of IFP
	:type tifp_type: str
	:return: table containing all the detected IFPs
	:rtype: pandas DataFrame
	:raises :class:'Exception': should there be a difference in the length of the IFPs stored in the same output file
	:raises :class:'ValueError': should be the type of IFP not implemented
	'''
	values = []
	header = None
	counter = 0

	if ifp_type == 'regular':
		interactions = ['HYD', 'FTF', 'ETF', 'HBD', 'HBA', 'CAT', 'ANI']
	elif ifp_type == 'polar':
		interactions = ['HBD', 'HBA', 'CAT', 'ANI', 'MCO']
	elif ifp_type == 'extended':
		interactions = ['HYD', 'FTF', 'ETF', 'HBD', 'HBA', 'CAT', 'ANI', 'PCI', 'MCO']
	else:
		raise Exception('Fingerprint format non-available')

	with open(file, 'r') as ifp_file:
		for line in ifp_file:
			if line.startswith('|') and not line.startswith('|ERROR'):
				if line.startswith('|WARNING'):
					values.append([np.nan for _ in range(len(values[0]))])
				elif header is None or header == line:
					header = line
				elif len(header) == len(line):
					print(f'Fingerprint header changed during execution at fingerprint {counter}, please check the output.\nThe fingerprint is still generated with the first detected fingerprint head')
				else:
					raise Exception(f'Inconsistent number of residues with previous fingerprints at line {counter}\nImpossible to create the output file')
				
			elif line.startswith('0') or line.startswith('1'):
				try:
					values.append([int(char) for char in line[:-1]])
				except ValueError as ve:
					pdb.set_trace()

	header_residues = header.split('|')
	header_residues[-1]=header_residues[-1][:-1]

	interaction_header = []
	

	for residue in header_residues[1:]:
		for code in interactions:
			interaction_header.append(f'{residue} {code}')

	if len(values[0]) != len(interaction_header):
		return pd.DataFrame(data = ['' for _ in range(len(interaction_header))])

	else:

		return pd.DataFrame(data = values, columns = interaction_header)
		

def filter_interaction(fingerprint, interactions):
	'''
	Filters the interaction type selecting only the bit associated to the selected interaction types.

	:param fingerprint: table containing IFPs
	:type fingerprint: pandas DataFrame
	:param interactions: interaction(s) to select from the IFP
	:type interactions: str of list of str
	:returns: filtered IFPs
	:rtype: pandas DataFrame
	:raises :class:'ValueError': Should the interaction not being between the recognized ones
	'''

	if not isinstance(interactions, list):
		interactions = [interactions]
    	
	recognized_interactions = ['HYD', 'FTF', 'ETF', 'HBD', 'HBA', 'CAT', 'ANI', 'PCI', 'MCO']
	inter_columns = list()

	for interaction in interactions:
		if interaction not in recognized_interactions:
			raise ValueError(f'Unrecognized interaction {interaction}')

	for index in fingerprint.columns:
		if index[-3:] in interactions:
			inter_columns.append(index)

	return fingerprint[inter_columns]




def filter_residues(fingerprint, residues):
	'''
	Filters the interaction type selecting only the bit associated to the selected residues.

	:param fingerprint: table containing IFPs
	:type fingerprint: pandas DataFrame
	:param residues: residue(s) to select from the IFP
	:type residues: str of list of str
	:returns: filtered IFPs
	:rtype: pandas DataFrame
	'''
	if not isinstance(residues, (list, np.ndarray)):
		residues = [residues]
	if not isinstance(fingerprint, pd.DataFrame):
		raise TypeError('The given fingerprint should be a DataFrame instance')
	columns = list()
	for bit in fingerprint.columns:
		bit_temp = bit.replace(' ', '')
		if int(bit_temp[2:-3]) in residues:
			columns.append(bit)

	return fingerprint[columns]