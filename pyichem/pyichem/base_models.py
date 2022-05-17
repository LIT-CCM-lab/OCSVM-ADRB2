import pandas as pd
import yaml
import os
import subprocess
import numpy as np



STRUCTURES_PATH = 'ichem_outputs/structures'

class IChem():

	'''
	The IChem class contains different type of calculation that can be performed by the IChem software suite.
	Functions of the class are used by all child classes.

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
	'''

	def __init__(self, receptor_mol2 = None, ligand_mol2 = None):
		'''Constructor method'''

		self.ligand_mol2=ligand_mol2
		self.receptor_mol2=receptor_mol2
		super().__init__()
		with open(f'{os.path.dirname(os.path.realpath(__file__))}/software_path.yml', 'r') as soft_p:
				self.software_paths = yaml.load(soft_p, Loader=yaml.FullLoader)
		self.ichem_path=self.software_paths['IChem']
		#self.ichem_path='IChem'

	def __check_files(self):

		'''
		Verifies that the correct output are given to the class before performing further calculations.
		It checkes the existance of the ligand and receptor files.
		It check that the same number of ligand and receptor file are given as input
		'''
		
		if self.receptor_mol2 is None:
			raise Exception("There are no mol2 files for the receptor")
		if self.ligand_mol2 is None:
			raise Exception("There are no mol2 files for the ligand")

		n_receptor = len(self.receptor_mol2)
		n_ligand = len(self.ligand_mol2)

		if not n_ligand == n_receptor:
			raise Exception("The same number of receptor and ligand files are necessary for the desired operation")

	def delete_mol2(self):
		'''Deletes all structure file present in the structure folder generated by the module'''
		subprocess.call(['rm', '-r', STRUCTURES_PATH])

	def preparation(self, path):
		'''
		Verifies that the given inputs are correct for the object.
		Generates the required folder to store the outputs

		:param path: folder to create to store the outputs
		:type path: str
		'''
		self.__check_files()

		if not os.path.isdir(path):
			os.makedirs(path)


class BatchCalculation(IChem):
	'''
	The BatchCalculation class is a child of the IChem class, it contains all the tools that can be called using an input file allowing for batch calculations.
	Batch calculations are significantly faster than running the single calculations.
	At the moment it contains only tools related to the detection and storing of interactions.

	:param receptor_mol2: files containing the structure of the protein to use in IChem calculations
	:type receptor_mol2: list of str, optional
	:param ligand_mol2: files containing the structure of the ligands to use in IChem calculations
	:type ligand_mol2: list of str, optional
	:param folder: folder where to store the output files
	:type folder: str
	:param output_p: prefix used to generate the names of the output files
	:type output_p: str
	:param software: indicate which tool to call from the IChem suite
	:type software: str
	:param output_s: suffix used to generate the names of the output files
	:type output_s: str, optional
	:param opt: it contains the various optional settings that can be given to the IChem tool
	:type opt: str, optional
	:param output_f: the tool generates an output file or not
	:type output_f: bool, optional
	:parm stdout_name: name of the file containg the stdout of IChem
	:type stdout_name: str, optional
	'''

	def __init__(self, receptor_mol2, ligand_mol2, folder, output_p, software, output_s = '', opt = '', output_f = True, stdout_name = 'ichem_stdout.txt'):
		'''Constructor method'''
		super().__init__(receptor_mol2, ligand_mol2)
		self.folder = folder
		self.stdout = stdout_name
		self.software = software
		self.output_prefix = output_p
		self.options = opt
		self.output_suffix = output_s
		self.output_file = output_f
		self.preparation(self.folder)

	def change_rules(self, parameters, new_values):
		'''
		Change the topological definition of protein ligand interactions from the default definitions in IChem.
		For a list of the code associated to each topological parameter consult IChem user manual.
		The function modifies the object attribute options.

		Accepted parameters:
		DHB, distance H-bond
		DHYD, distance hydrophobic contact
		DIO, dinstance ionic interaction
		DME,
		DAR, distance aromatic
		DPIC, dinstance pi-cation
		AH
		ATH
		AARFF
		ATARFF
		AAERF
		ATAREF
		APIC
		ATPIC


		:param parameters: topological parameters to change
		:type parameters: str or list of str
		:param new_values: new values for the indicated topological parameters
		:type new_values: float or list of floats
		'''
		rule_commands = {'DHB' : '–D_Hb', 'DHYD' : '–D_Hyd', 'DIO': '–D_Io', 'DME': '–D_Me', 'DAR': '-D_Ar',
		 'DPIC': '-D_PIC', 'AH': '-a_H', 'ATH': '-at_H', 'AARFF' : '–a_ArFF', 'ATARFF' : '–at_ArFF',
		  'AAREF' : '–a_ArEF', 'ATAREF' : '–at_AREF', 'APIC': '–a_Pic', 'ATPIC' : '–at_PIC' }

		accepted_parameters = list(rule_commands)

		new_rule = ''

		if not isinstance(parameters, list):
			parameters = [parameters]
		if not isinstance(new_values, list):
			new_values = [new_values]

		if len(parameters) != len(new_values):
			raise ValueError('Inconsistent size between the parameters and the new values')

		for i, new_p in enumerate(parameters):
			if new_p in accepted_parameters:
				new_rule = new_rule + f'{rule_commands[new_p]} {new_values[i]} '
			else:
				raise ValueError('Uknown parameter set for the change of rules')

		self.options = new_rule + self.options

	def _write_input(self, input_file):
		'''
		Write the input file given to IChem.

		:param input_file: name of the input file
		:type input_file: str
		'''

		n_files = len(os.listdir(self.folder))
		output_location = []
		with open(f'{self.folder}/{input_file}', 'w') as input_f:
			for i, (receptor, ligand) in enumerate(zip(self.receptor_mol2, self.ligand_mol2)):
				file_number = int((i+n_files))
				if self.output_file:
					input_f.write(f'{self.options} {self.software} {receptor} {ligand} {self.output_prefix}{str(file_number)}{self.output_suffix}\n')
					output_location.append(f'{self.output_prefix}{str(file_number)}{self.output_suffix}')

				else:
					input_f.write(f'{self.options} {self.software} {receptor} {ligand} \n')


		if len(output_location) != 0:
			self.output_location = output_location
		self.input_file = input_file

	def calculate(self, input_file = 'ichem_input.in', stdout_c = True, stderr_c = True):
		'''
		Generates the input file, launches the calculation, and saves stdout and stderr to files.

		:param input_file: name of the input file
		:type input_file: str
		:param stdout_c: capture the stdout stream from IChem
		:type stdout_c: bool
		:param stderr_c: capture the stderr stream from IChem
		:type stderr_c: bool
		'''

		
		self._write_input(input_file)

		try:
			process_output = subprocess.run([self.ichem_path, '-F', f'{self.folder}/{self.input_file}'],
			 	stdout = subprocess.PIPE if stdout_c else None,
			 	stderr = subprocess.PIPE if stderr_c else None)
		except FileNotFoundError:
			raise FileNotFoundError(f'The selected IChem path was incorrect: {self.ichem_path}.\nPlease update the file at: {os.path.dirname(os.path.realpath(__file__))}/software_path.yml')

		if process_output.stdout is not None:
			with open(self.stdout, 'w') as out:
				out.write(process_output.stdout.decode('utf-8'))

		if process_output.stderr is not None:
			with open('ichem_stderr.txt', 'a') as err:
				err.write(process_output.stderr.decode('utf-8'))

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
			for line in inpt:
				command = line.split()
				command.insert(0, self.ichem_path)
				try:
					process_output = subprocess.run(command,
			 			stdout = subprocess.PIPE if stdout_c else None,
			 			stderr = subprocess.PIPE if stderr_c else None)
				except FileNotFoundError:
					raise FileNotFoundError(f'The selected IChem path was incorrect: {self.ichem_path}.\nPlease update the file at: {os.path.dirname(os.path.realpath(__file__))}/software_path.yml')

				if process_output.stdout is not None:
					with open(self.stdout, 'a') as out:
						out.write(process_output.stdout.decode('utf-8'))

				if process_output.stderr is not None:
					with open('ichem_stderr.txt', 'a') as err:
						err.write(process_output.stderr.decode('utf-8'))

	def map_results(self, name = False, energy = False):
		'''
		Generates a pandas DataFrame containing all output files generated by IChem and the protein and ligand files used to generate it.
		It has the option of reading the mol2 files looking for the name of the molecule and its energy
		:param name: Search for the name of the molecule in the mol2 file
		:type name: bool
		:param energy: Search for the docking score in the mol2 file
		:type energy: bool
		:returns: table containing the file used to generate the output file and the output file, if selected also returns the enrgy and docking score
		:rtype: pandas DataFrame
		'''

		return pd.DataFrame(data = np.array([self.receptor_mol2, self.ligand_mol2, self.output_location]).T, columns = ['Receptor_file', 'Ligand_file', 'Output_file'])

		

	def read_input_file(self, input_file = 'ichem_input.in'):
		'''
		Reads the previously generated input file to obtain information about the input and output files

		:param input_file: file containing the inputs generated by IChem
		:type input_file: str, optional
		'''
		with open(f'{self.folder}/{input_file}') as inp:
			self.receptor_mol2 = list()
			self.ligand_mol2 = list()
			self.output_location = list()
			for line in inp:
				elements = line.split()
				for i, part in enumerate(elements):
					if part == self.software:
						starting_index = i+1
						break
				self.receptor_mol2.append(elements[starting_index])
				self.ligand_mol2.append(elements[starting_index+1])
				self.output_location.append(elements[starting_index+2])

	def read_map_file(self, map_file):
		map_df = pd.read_csv(map_file)

		self.receptor_mol2 = map_df['Receptor_file'].tolist()
		self.ligand_mol2 = map_df['Ligand_file'].tolist()
		self.output_location = map_df['Output_file'].tolist()