'''
Created by Patrick Kilecdi

Certain utilities called by the other functions in section
'''

def aims_extract_energy(path):
	'''
	Extracts energy from a given aims output file path
	'''
	aims_out = open(path)
	while True:
		line = aims_out.readline()
		if not line: return False  # energy not converged
		if '| Total energy of the DFT /' in line:
			tokens = line.split()
			energy = float(tokens[11])  # converts from SI string to float
			return energy

