from merlin import *

import matplotlib.pyplot as plt
import numpy as np

def test_functional_form():
	eps = 1.0
	sig = 3.0
	lj = lambda r: 4*eps*( (sig/r)**12 - (sig/r)**6 )
	r_cutoff = 4.
	cut = lambda r: (r**6 + r_cutoff**6)**(1./6)

	lj2 = lambda r: 4*eps*( (sig/cut(r))**12 - (sig/cut(r))**6 )

	xx = np.arange(2.,10.,0.1)
	Es = [lj(r) for r in xx]
	E2s = [lj2(r) for r in xx]

	plt.plot(xx, Es)
	plt.plot(xx, E2s)
	plt.ylabel('E')
	plt.ylim(-10.0, 10.0)
	plt.show()

def read_old_gaussian():
	source_directory = '/fs/home/wmc62/Documents/perovskites/smrff/gaussian'
	for root, dirs, file_list in os.walk(source_directory):
		for name in file_list:
			if name.endswith('.log'):
				if not (name.startswith('PbCl2') or name.startswith('PbCl_24')): continue
				if not '_vac' in name: continue
				if '_new_' in name : continue
				if 'bad' in name : continue
				if name.startswith('PbCl2t'): continue
				result = g09.parse_atoms(source_directory+'/'+name, check_convergence=True)
				if result:
					energy, atoms = result
					name = name[:-4]
					print name, '\t', ' '.join([a.element for a in atoms])
					if os.path.isdir('orca/'+name):
						#orca.job(name, '! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid3 FinalGrid5 SlowConv', atoms, queue=None, grad=True, previous=name).wait()
						for a in atoms:
							if a.element=='Pb': a.label='907'
							if a.element=='Cl': a.label='344'
						files.write_cml(atoms, name='orca/'+name+'/system.cml')
read_old_gaussian()

