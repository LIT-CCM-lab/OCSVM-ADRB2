import argparse
import sys
import pdb
import pytraj as pt
import subprocess
import os
import time
import numpy as np
import mol2_trajectory
from random import randint

def print_progress (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r", refreshRate=0.01):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()




def align_traj(traj_inp, mask_align, ref, mask_ref):

	#pdb.set_trace()

	memory_traj=pt.Trajectory.from_iterable(traj_inp)
	memory_traj=memory_traj.autoimage()
	memory_traj=memory_traj.superpose(mask=f':{mask_align}', ref=ref, ref_mask=f':{mask_ref}')

	
	
	rmsd = pt.rmsd(memory_traj, mask=f':{mask_align}@CA', ref=ref, ref_mask=f':{mask_ref}@CA', nofit = True)
	max_rmsd = np.amax(rmsd)
	if max_rmsd > 10:
		memory_traj = None

	return memory_traj


def main(args):

	#pdb.set_trace()
	report = list()
	report.append(f'Trajectory to mol2 conversion report\n')
	report.append(f'Trajectory files: {args.trajectory[0]}')
	if len(args.trajectory) > 1:
		for tr in args.trajectory[1:]:
			report.append(f'\t{tr}')
	report.append(f'Topology file: {args.topology}')
	
	print('Loading trajectory ...')

	if args.last_frame is not None and args.skip_frames is not None:
		raise Exception('It is not possible to perform both stride and slicing on the same trajectory')
	elif args.last_frame is not None:
		traj = pt.iterload(args.trajectory, args.topology, frame_slice = (args.first_frame, args.last_frame))
		report.append(f'Starting trajectory frame: {args.first_frame}\nLast trajectory frame: {args.last_frame}')
	else:
		traj = pt.iterload(args.trajectory, args.topology, stride = args.skip_frames)
		report.append(f'Trajectory stride: {args.skip_frames}')

	reference = pt.iterload(args.reference)

	print(f'Trajectory loaded with {traj.n_frames} frames')
	report.append(f'Loaded frames: {traj.n_frames}')

	stability_frame = None
	memory_traj = None
	error_message = ''

	if args.alignment_ref is None:
		args.alignment_ref = args.alignment

	if traj.n_frames < args.chunk:
		chunk = traj.n_frames
	else:
		chunk = args.chunk

	for i, traj_i in enumerate(traj.iterchunk(chunksize=chunk)) :
		counter = 0
		print_progress(i*chunk, traj.n_frames)
		
		while memory_traj is None:
			memory_traj=align_traj(traj_i, args.alignment, reference, args.alignment_ref)
			counter +=1
			if counter > 100:
				memory_traj=align_traj(traj_i, args.alignment, reference, args.alignment_ref, None)
				error_message += f'The program was not able to correctly align frames {i*chunk} to {(i+1)*chunk-1}\n'
				break
		
		 
		mol2_traj=mol2_trajectory.Trajectory(traj = traj_i, receptor_mask = args.receptor, ligand_mask = args.ligand, ff = args.force_field, pdb = args.pdb_conversion)
		mol2_traj.write_mol2()


		stability_frame = memory_traj[-1]
		memory_traj = None

		

	print_progress(traj.n_frames, traj.n_frames)

	print('Completing conversion to mol2 files ...')
	mol2_traj=mol2_trajectory.Trajectory(ff = args.force_field, pdb = args.pdb_conversion)
	mol2_traj.load_mol2()
	print('Fixing receptor mol2 files\n')
	mol2_traj.fix_receptor_file()
	print('\n')
	print('Fixing ligand mol2 files\n')
	mol2_traj.fix_ligand_file()
	print('\n')

	if error_message == '':
		print('Process completed without errors')
	else:
		print('Process completed')
		print(error_message)
		report.append(error_message)

	with open(args.report, 'w') as rep:
		rep.writelines('\n'.join(report))


		
	

if __name__ == "__main__":
	parser=argparse.ArgumentParser()
	parser.add_argument('-top', '--topology', help='input topology file', required = True)
	parser.add_argument('-traj', '--trajectory', nargs='+', help='input trajectory file(s)', required = True)
	parser.add_argument('-r', '--receptor', help='Residues of the receptor', required = True)
	parser.add_argument('-l', '--ligand', help='Residues of the ligand', required = True)
	parser.add_argument('-ref', '--reference', help='Reference file', required = True)
	parser.add_argument('-a', '--alignment', help='Residues to use for the alignment', required = True)
	parser.add_argument('-ar', '--alignment_ref', default=None, help='Residues to use for the alignment of the reference structure')
	parser.add_argument('-sf', '--skip_frames',default=None, help='Skip this many frames')
	parser.add_argument('-p', '--pdb_conversion',default=None, help='.csv file containing the atom types')
	parser.add_argument('-if', '--first_frame', default = 0, help = 'frame from which start the conversion to mol2 file', type = int)
	parser.add_argument('-lf', '--last_frame', default = None, help = 'frame where the conversion to mol2 file ends', type = int)
	parser.add_argument('-ff', '--force_field', default = None, help = 'force field in which the trajectory atom types are defined.\n Set to charmm if CHARMM is used, default option considers AMBER atom types')
	parser.add_argument('-c', '--chunk', default = 100, help = 'Number of frames to be processed at each iteration of the converter')
	parser.add_argument('-rep', '--report', default = 'trajectory_conversion_report.txt', help = 'Name of the trajectory conversion file')
	

	parser.set_defaults(func=main)
	args=parser.parse_args()
	status = args.func(args)
	sys.exit(status)