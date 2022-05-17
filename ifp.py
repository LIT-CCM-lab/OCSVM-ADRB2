import argparse
import sys
import pdb
import subprocess
import os
import numpy as np
from mol2_trajectory import Trajectory
from pyichem import ifp

def main(args):
	
	mol2_traj=Trajectory()
	mol2_traj.load_mol2()

	print('Detecting protein-ligand interactions')

	fingerprint = ifp.Ifp(mol2_traj.receptor_mol2, mol2_traj.ligand_mol2)
	fingerprint.change_rules(['DAR'], [5.0])
	fingerprint.calculate_lbl()

	fingerprint.read_ifp()
	print('Saving IFP to csv file\n')
	fingerprint.fingerprints.to_csv('ifp.csv', index = False)
	fingerprint.map_results().to_csv('ifp_map.csv', index = False)
	print('IFP saved')
	

		
	

if __name__ == "__main__":
	parser=argparse.ArgumentParser()

	parser.set_defaults(func=main)
	args=parser.parse_args()
	status = args.func(args)
	sys.exit(status)