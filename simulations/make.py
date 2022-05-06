import os, sys
import re
sys.path.append('/home/nam4')
import numpy as np
import orbifold_tools.rigid.lammps as rs

def make(
	output,
	init_config, 
	pot_filename,
	cpu_count, 
	rec_freq, 
	n_atom_types, 
	T, 
	np_width,
	n_equil,
	n_prod,
	Tdamp=1, 
	dt=0.005, 
	bins=10000
	):

	raw = open('/home/nam4/rigid_np/run.lmp.template', 'r').read()
	data = re.sub('__INIT_CONFIG__', str(init_config), raw)
	data = re.sub('__CPU_COUNT__', str(cpu_count), data)	
	data = re.sub('__REC_FREQ__', str(rec_freq), data)
	data = re.sub('__RES_FREQ__', str(int(rec_freq*100)), data)
	data = re.sub('__NATOM__', str(n_atom_types), data)
	data = re.sub('__TEMP__', str(T), data)
	data = re.sub('__TDAMP__', str(Tdamp), data)
	data = re.sub('__TIMESTEP__', str(dt), data)
	data = re.sub('__TABLE_BINS__', str(bins), data)
	data = re.sub('__NP_WIDTH__', str(np_width), data)
	data = re.sub('__N_EQUIL__', str(int(n_equil)), data)
	data = re.sub('__N_PROD__', str(int(n_prod)), data)

	command = ""
	for i in range(1, n_atom_types + 1):
		for j in range(i, n_atom_types + 1):
			if i == j and i < n_atom_types:
				command += "pair_coeff\t{} {} {} {}_{}\n".format(i, j, pot_filename, 1, 2)
			else:
				command += "pair_coeff\t{} {} {} {}_{}\n".format(i, j, pot_filename, i, j)
	data = re.sub('__PAIR_COEFF_WCA__', command, data)

	command = ""
	for i in range(1, n_atom_types + 1):
		for j in range(i, n_atom_types + 1):
			command += "pair_coeff\t{} {} {} {}_{}\n".format(i, j, pot_filename, i, j)
	data = re.sub('__PAIR_COEFF__', command, data)

	open(output, 'w').write(data)

if __name__ == '__main__':
	import sys
	coords, ids, np_width, core_sigma, sigma = rs.extract_xyz(sys.argv[1])
	n1, n2 = int(sys.argv[2]), int(sys.argv[3]) # Number of each chirality
	rho = float(sys.argv[4]) # System total coverage / packing fraction
	Lx = Ly = np.sqrt(np.pi*np_width**2*(n1+n2)/rho)
	print('np_width = {}'.format(np_width))

	if len(sys.argv) > 7:
		buffer_ = float(sys.argv[7])
	else:
		buffer_ = sigma

	init_filename = 'init.lammps'
	rs.create_init_config(coords, ids, box=[Lx,Ly], n=[n1,n2], sigma=buffer_,
		filename=init_filename)

	pot_filename = 'pair.potentials'
	table_bins = 10000
	rs.tabulate_potentials(pot_filename,
                        table_bins=table_bins,
                        n_types=len(np.unique(ids)),
                        sigma=sigma,
                        core_sigma=core_sigma,
                        alpha=9,
                        r_cut=3.0)

	n_atom_types = len(np.unique(ids))

	make(
		output='run.lmp',
		init_config=init_filename,
		pot_filename=pot_filename,
		cpu_count=int(sys.argv[6]), 
		rec_freq=10000, 
		n_atom_types=n_atom_types, 
        	T=float(sys.argv[5]), 
		np_width=np_width,
		n_equil=1e6,
		n_prod=50e6,
		Tdamp=10, 
		dt=0.005, 
		bins=table_bins
	)
