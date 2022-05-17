import argparse
import sys
import pdb
import subprocess
import os
import numpy as np
from mol2_trajectory import Trajectory
from pyichem import ints

def main(args):
	mol2_traj=Trajectory()
	mol2_traj.load_mol2()

	print('Starting calculations')

	if args.default:
		print('Calculating interactions using the default definitions')
		fingerprint = ints.Ints(mol2_traj.receptor_mol2, mol2_traj.ligand_mol2, type_int = 'MERG', new_hyd = True)
		fingerprint.change_rules(['DAR'], [5.0])
		fingerprint.calculate(input_file = 'interactions_newhyd.in')
		fingerprint.map_results().to_csv('interactions_map_newhyd.csv', index = False)

	if args.new:
		print('Calculating interactions using the Newhyd defintion of hydrophobic contacts')
		fingerprint = ints.Ints(mol2_traj.receptor_mol2, mol2_traj.ligand_mol2, type_int = 'MERG', new_hyd = False)
		fingerprint.change_rules(['DAR'], [5.0])
		fingerprint.calculate(input_file = 'interactions.in')
		fingerprint.map_results().to_csv('interactions_map.csv', index = False)

	print('Calculation completed')
	

		
	

if __name__ == "__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('-d', '--default', default = True, help = 'Compute IPA with default defintion of hydrophobic contacts')
	parser.add_argument('-n', '--new', default = True, help = 'Compute IPA with the Newhyd defintion of hydrophobic contacts')

	parser.set_defaults(func=main)
	args=parser.parse_args()
	status = args.func(args)
	sys.exit(status)