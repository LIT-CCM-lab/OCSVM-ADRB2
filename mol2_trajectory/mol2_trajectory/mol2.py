
substructure_card = '@<TRIPOS>SUBSTRUCTURE\n'
atom_card= '@<TRIPOS>ATOM\n'
bond_card = '@<TRIPOS>BOND\n'
comment_card = '@<TRIPOS>COMMENT\n'
molecule_card = '@<TRIPOS>MOLECULE\n'

separators = [molecule_card, atom_card, bond_card, substructure_card, comment_card]


class mol2_file():
	def __init__(self,in_file, pdb_conversion = None, backbone_tag = False, ff_conversion = None, out_file = None, smiles = None):
		self.file = in_file
		self.check_pdb_conversion = pdb_conversion
		self.backbone_tag = backbone_tag
		self.ff_conversion = ff_conversion
		self.pdb_conversion = pdb_conversion
		self.smiles = None

		if out_file is None:
			self.out_file = in_file
		else:
			self.out_file = out_file


	def fix_mol2(self):
		self.generate_blocks()
		self.fix_atom_block()
		self.fix_molecule_block()
		self.fix_bond_block()
		self.fix_substructure_block()
		self.write_mol2()

	


	def generate_blocks(self):

		with open(self.file, 'r') as output:
			data = output.read()

		parts = ['',data]

		blocks = list()

		for sep in separators:
			parts = parts[-1].split(sep)
			blocks.append(parts[0])

		if len(parts) == 2:
			blocks.append(parts[1])

		del blocks[0]

		self.blocks = blocks

	def write_mol2(self):
		
		with open(self.out_file, 'w') as output:
			for card, block in zip(separators, self.blocks):
				output.writelines(card)
				for line in block:
					output.writelines(line)
				

	def fix_molecule_block(self):
		lines = self.blocks[0].split('\n')
		lines[0] = 'mol2 file generated by mol2_trajectory'
		lines = [line+'\n' for line in lines]
		self.blocks[0] = lines

	def fix_atom_block(self):
		lines = self.blocks[1].split('\n')
		new_lines = list()

		self.aromatics = list()
		self.sp2 = list()
		self.backbone = list()
		self.amide = list()

		
		if self.smiles is None:
			for line in lines:
				if len(line) != 0:
					new_line, _ = self.fix_atoms_line(line)
					new_lines.append(new_line)
		else:
			smiles_ref = 0
			for line in lines:
				if len(line) != 0:
					new_line, smiles_ref = self.fix_atoms_line(line, smiles_ref = smiles_ref)
					new_lines.append(new_line)
				

		self.blocks[1] = new_lines

	def fix_bond_block(self):
		lines = self.blocks[2].split('\n')
		new_lines = list()

		for line in lines:
			new_lines.append(self.fix_bond_line(line))

		self.blocks[2] = new_lines

	def fix_substructure_block(self):
		lines = self.blocks[3].split('\n')
		new_lines = list()

		for line in lines:
			new_lines.append(self.fix_substructure_line(line))

		self.blocks[3] = new_lines


	def fix_atoms_line(self, line, smiles_ref = None):
		'''
		Standardize mol2 atom block for use with IChem

		:param line: line from the mol2 file
		:type line: str
		:param conversion_file: conversion between a force field and SYBYL atom types
		:type conversion-file: dict

		:return: the standardized line, atom number if the atom is aromatic, atom number if the atom is an O or N sp2
		:rtype: str, str, str
		'''
		
		sybyl_aro = ['C.ar', 'N.ar', 'O.co2', 'N.pl3', 'C.cat'] #O.co2, N.pl3, C.cat are considered as aromatic in order to implement the ar bond type.
		sybyl_sp2 = ['O.2', 'N.2', 'C.2', 'S.2']
		backbone = ['C', 'CA', 'O', 'N', 'H', 'HA', 'HA1', 'HA2', 'HA3']
		amide_atoms = ['C', 'N']

		trp_atoms = {'CG': 'C.ar'}
		hid_atoms = {'ND1': 'N.pl3', 'NE2': 'N.ar'}
		hie_atoms = {'ND1': 'N.ar', 'NE2': 'N.pl3'}
		hip_atoms = {'ND1': 'N.pl3', 'NE2': 'N.pl3'}
		arg_atoms = {'NE': 'N.pl3', 'NH1': 'N.pl3', 'NH2': 'N.pl3', 'CZ': 'C.cat'}

		line = line.replace('S.O', 'S.o')
		line = line.replace('S.O2', 'S.o2')
		components = line.split()

		
		if len(components) == 0:
			return line
		if len(components) < 9:
			#pdb.set_trace()
			components.insert(1, components[0][-4:])
			components[0] = components[0][:-4]

			
		if self.pdb_conversion is not None:
			components[5] = to_atom_type(components[1], components[7], self.pdb_conversion)
		if self.ff_conversion is not None:
			components[5] = to_sybyl(components[5], self.ff_conversion)

		if len(components) == 10:
			status_bit = components[9]
		elif components[1] in backbone and self.backbone_tag:
			status_bit = 'BACKBONE'
			self.backbone.append(components[0])
				
			if components[1] in amide_atoms:
				self.amide.append(components[0])
		else:
			status_bit = ''
		'''
		if smiles_ref is not None and components[5] != 'H':
			#pdb.set_trace()
			if not components[5].upper().startswith(self.smiles[smiles_ref].upper()):
				print(components[5])
				print(self.smiles)
				raise Exception('The SMILES order of atoms and the mol2 file order of atoms are not compatible')
			elif self.smiles[smiles_ref] == 'c':
				components[5] = 'C.ar'
			elif self.smiles[smiles_ref] == 'n':
				components[5] = 'N.ar'

			smiles_ref += 1
		'''

		if components[7][:3] == 'TRP':
			if components[1] in list(trp_atoms.keys()):
				components[5] = trp_atoms[components[1]]

		elif components[7][:3] == 'HID' or components[7][:3] == 'HSD':
			if components[1] in list(hid_atoms.keys()):
				components[5] = hid_atoms[components[1]]

		elif components[7][:3] == 'HIE' or components[7][:3] == 'HSE':
			if components[1] in list(hie_atoms.keys()):
				components[5] = hie_atoms[components[1]]

		elif components[7][:3] == 'HIP':
			if components[1] in list(hip_atoms.keys()):
				components[5] = hip_atoms[components[1]]

		elif components[7][:3] == 'ARG':
			if components[1] in list(arg_atoms.keys()):
				components[5] = arg_atoms[components[1]]

		elif components[7][:3] == 'HIS':
			components[5] = 'N.ar' if components[5] == 'N.2' else components[5]
			components[5] = 'N.pl3' if components[5] == 'N.3' else components[5]



		components[7] = check_res(components[7][:3])+components[7][3:]




		if len(components[7]) == 3:
			components[7] = components[7]+components[6]


		new_line = f'{components[0]:>7} {components[1]:>5}{components[2]:>14}{components[3]:>10}{components[4]:>10} {components[5]:<11}{components[6]} {components[7]:<8}{components[8]:>9} {status_bit}\n'

		if components[5] in sybyl_aro:
			self.aromatics.append(components[0])
		elif components[5] in sybyl_sp2:
			self.sp2.append(components[0])

		return new_line, smiles_ref

	def fix_bond_line(self, line):
		'''
		Standardize mol2 bond block for use with IChem.
		Converts single bond to correct aromatic or double bonds

		:param line: line from the mol2 file
		:type line: str
		:param aromatics: aromatic atoms
		:type aromatics: list of str
		:param sp2: sp2 atoms
		:type sp2: list of str

		:return: the standardized line
		:rtype: str
		'''
		components = line.split()

		backbone_bit = ''

		if len(components) == 0:
			return line

		if components[1] in self.backbone and components[2] in self.backbone:
			if components[1] in self.amide and components[2] in self.amide:
				components[3] = 'am'
				backbone_bit = '  BACKBONE|INTERRES'
			else:
				backbone_bit = '  BACKBONE'

		elif components[1] in self.aromatics and components[2] in self.aromatics:
			components[3] = 'ar'
		elif components[1] in self.sp2 and components[2] in self.sp2:
			components[3] = 2
		elif (components[1] in self.aromatics or components[1] in self.sp2) and (components[2] in self.aromatics or components[2] in self.sp2):
			components[3] = 'ar'


		return f'{components[0]:>6}{components[1]:>5}{components[2]:>5} {components[3]}{backbone_bit}\n'

	def fix_substructure_line(self, line):
		'''
		Standardize mol2 substructure block for use with IChem

		:param line: line from the mol2 file
		:type line: str

		:return: the standardized line
		:rtype: str
		'''
		components = line.split()
		if len(components) >= 7:
			new_line = f'{components[0]:>7}{check_res(components[1][:3])+components[0]:>7}{components[2]:>16} RESIDUE{4:>14} A {check_res(components[1][:3])}{0:>6}\n'
			return new_line
		else:
			return line+'\n'

def to_sybyl(atom, conversion_dict):
	'''
	Converts the atom type from a force field to SYBYL.

	:param atom: ff atom type
	:type atom: str
	:param conversion_dict: conversion between ff and SYBYL
	:type conversion_dict: dict

	:return: Converted atom type
	:rtype: str

	:raises: :class:'Exception': Should the ff atom type not being in the conversion file
	'''
	for sybyl in list(conversion_dict):
		if atom in conversion_dict[sybyl]:
			return sybyl
	raise Exception(f'atom type {atom} is not supported, please update the conversion file.')

def to_atom_type(atom_name, res_name, conversion_file):
	'''
	Finds the atom type of a specific atom from a list.

	:param atom_name: name of the atom
	:type atom_name: str
	:param res_name: name of the residue
	:type res_name: str
	:param conversion_dict: conversion between atom and residue and the atom type
	:type conversion_dict: pandas DataFrame

	:return: Converted atom type
	:rtype: str

	'''
	
	return(conversion_file[res_name][atom_name])

def check_res(res):
	'''
	Converts non-standard residue names into standard residue names

	:param res: residue name
	:type res: str
	:returns: standard residue name
	:rtype: str
	'''
	res_name_dict = {'CYX': 'CYS', 'HIE': 'HIS', 'HIP': 'HIS', 'HID': 'HIS', 'HSD': 'HIS', 'HSE': 'HIS', 'ASH': 'ASP', 'LYN': 'LYS', 'GLH': 'GLU'}

	if res in res_name_dict:
		return res_name_dict[res]
	else:
		return res