import os.path
import yaml

path = os.path.dirname(os.path.realpath(__file__)) + '/conversion_file/'

def load_ff(ff = 'charmm'):
	if ff == 'charmm':
		with open(f'{path}charmm2sybyl.yml') as c2s:
			return yaml.load(c2s, Loader=yaml.FullLoader)
	elif ff == 'amber':
		with open(f'{path}amber2sybyl.yml') as c2s:
			return yaml.load(c2s, Loader=yaml.FullLoader)
	elif ff is None:
		return None
	else:
		raise ValueError(f'Unrecognized forcefield {ff}')

	return ff_conversion

def load_pdb_c(pdb_c):
		if pdb_c is not None:
			with open(pdb_c) as p:
				return yaml.load(p, Loader=yaml.FullLoader)
		else:
			return None

def get_path_files(folder):
		if os.path.isdir(folder):

			list_file=[]

			files = os.listdir(folder)

			for f in files:
				list_file.append(f'{folder}/{f}')
			return sorted(list_file)

		else:
			return None

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