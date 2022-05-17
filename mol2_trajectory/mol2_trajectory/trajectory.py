import numpy as np
import os
import pandas as pd
import pytraj as pt
import sys
sys.path.append('/projects/cxcr4/cxcr4/git_scripts/pyChem')
from mol2_trajectory.utils import load_ff, load_pdb_c, get_path_files, print_progress
from mol2_trajectory.mol2 import mol2_file

LIGAND_PATH='ichem_outputs/structures/ligand'
RECEPTOR_PATH = 'ichem_outputs/structures/receptor'
STRUCTURES_PATH = 'ichem_outputs/structures'



class Trajectory():
	'''
	Class containing a trajectory object.
	The class converts the trajectory into multiple .mol2 files to be used with IChem.
	The class supports amber trajectories by default and was adapted to handle CHARMM trajectories and Amber trajectories with a different conversion.

	:param traj: trajectory
	:type traj: pytraj Trajectory, optional
	:param receptor_mask: residues forming the protein
	:type receptor_mask: str, optional
	:param ligand_mask: residues forming the ligand
	:type ligand_mask: str, optional
	:param ff: the force field used for the trajectory
	:type ff: string, optional
	'''
	def __init__(self, traj = None, receptor_mask = None, ligand_mask = None, ff = None, pdb = None, ref_smiles = None):
		'''Constructor method'''
		self.traj = traj
		self.receptor_mask = receptor_mask
		self.ligand_mask = ligand_mask
		self.ff = ff
		self.pdb = pdb
		if isinstance(ref_smiles, list):
			self.ref_smiles = ref_smiles
		elif isinstance(ref_smiles, str):
			self.ref_smiles = [ref_smiles]
		elif ref_smiles is None:
			self.ref_smiles = None
		else:
			raise ValueError('The reference SMILES must be a string or list of strings')

	def write_mol2(self):
		'''Writes the trajectory to mol2 files'''

		if not os.path.isdir(LIGAND_PATH):
			os.makedirs(LIGAND_PATH)
		if not os.path.isdir(RECEPTOR_PATH):
			os.makedirs(RECEPTOR_PATH)
				
		if self.receptor_mask is not None:
			self.receptor_mol2 = self.__mol2_writer(self.receptor_mask, RECEPTOR_PATH, 'receptor_output')
	
		if self.ligand_mask is not None: 
			self.ligand_mol2 = self.__mol2_writer(self.ligand_mask, LIGAND_PATH, 'ligand_output')

	def load_mol2(self, receptor_folder = RECEPTOR_PATH, ligand_folder = LIGAND_PATH):
		'''
		Loads mol2 files into the trajectory

		:param receptor_folder: folder in which the receptor mol2 files are stored
		:type receptor_folder: str, optional
		:param ligand_folder: folder in which the ligand mol2 files are stored
		:type ligand_folder: str, optional
		'''
		self.receptor_mol2 = get_path_files(receptor_folder)
		self.ligand_mol2 = get_path_files(ligand_folder)

	def fix_receptor_file(self):
		'''The receptor files are changed in order to conform to the standards required by IChem'''

		self.__fix_file(self.receptor_mol2, True)

	def fix_ligand_file(self, protein = False):
		'''The ligand files are changed in order to conform to the standards required by IChem'''

		self.__fix_file(self.ligand_mol2, protein, smiles = self.ref_smiles)


	def __mol2_writer(self, mask, folder, name):
		'''
		Writes the trajectory to mol2 files.
		The files are renamed and standardised to be readable by IChem.

		:param mask: residues to convert as a mol2 file
		:type mask: str
		:param folder: folder where to store the output diles
		:type folder: str
		:param name: name of the output file
		:type name: str

		:returns: files containing the converted trajectory
		:rtype: list of str
		'''
		n_files = len(os.listdir(folder))
		target=pt.strip(self.traj, f"!(:{mask})")
		target_files = target.n_frames

		if self.ff is not None:
			pt.write_traj( f"{folder}/{name}.mol2",traj=target, 
			format='mol2', options='multi', overwrite=True)
		else:
			pt.write_traj( f"{folder}/{name}.mol2",traj=target, 
				format='mol2', options='multi sybyltype', overwrite=True)
		
		for i in range(1, target_files+1):
			os.rename(f"{folder}/{name}.mol2.{str(i)}", 
				f"{folder}/{name}_{str(i+n_files)}.mol2")
			

		return get_path_files(folder)

	def __fix_file(self, files, backbone, smiles = None):
		'''
		Standaridize the mol2 file to use it with IChem

		:param file: file to standardize
		:type file: str
		'''
		ff_conversion = load_ff(self.ff)
		pdb_conversion = load_pdb_c(self.pdb)

		if smiles is None:
			for i,file in enumerate(files):
				print_progress(i,len(files))
				to_fix_file = mol2_file(file, ff_conversion = ff_conversion, pdb_conversion = pdb_conversion, backbone_tag = backbone)
				to_fix_file.fix_mol2()
		else:
			if len(smiles) == 1:
				for file in files:
					to_fix_file = mol2_file(file, ff_conversion = ff_conversion, pdb_c = pdb_conversion, backbone_tag = backbone, smiles = smiles[0])
					to_fix_file.fix_mol2()
			else:
				for file, smiles_f in zip(files, smiles):
					to_fix_file = mol2_file(file, ff_conversion = ff_conversion, pdb_c = pdb_conversion, backbone_tag = backbone, smiles = smiles_f)
					to_fix_file.fix_mol2()