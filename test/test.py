from merlin import *

import matplotlib.pyplot as plt
import numpy as np
import time

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
				print name
				continue
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
					if not os.path.isdir('orca/'+name):
						if len(atoms)>3: continue
						orca.job(name, '! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid3 FinalGrid5 SlowConv', atoms, queue=None, grad=True).wait()
						for a in atoms:
							if a.element=='Pb': a.label='907'
							if a.element=='Cl': a.label='344'
						files.write_cml(atoms, name='orca/'+name+'/system.cml')

def read_old_gaussian_pairs():
	source_directory = '/fs/home/wmc62/Documents/perovskites/smrff/gaussian'
	for root, dirs, file_list in os.walk(source_directory):
		for name in file_list:
			if name.endswith('.log'):
				if not name.startswith('PbI+'): continue
				if name.endswith('SVP.log'): continue
				if '_' not in name: continue
				result = g09.parse_atoms(source_directory+'/'+name, check_convergence=True)
				if result:
					energy, atoms = result
					name = name[:-4]
					name = name.replace('I','Cl')
					for a in atoms: 
						if a.element=='I': a.element='Cl'
					print name, '\t', ' '.join([a.element for a in atoms])

					orca.job(name, '! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid3 FinalGrid5 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='1 1').wait()
					for a in atoms:
						if a.element=='Pb': a.label='907'
						if a.element=='Cl': a.label='344'
					files.write_cml(atoms, name='orca/'+name+'/system.cml')

def maintain_queue_size(queue, N):
	while len(queue)>=N:
		for job in queue:
			if job.poll() != None: #job is done
				queue.remove(job)
		time.sleep(0.1)

def pairwise_pbcl():
	x = 2.0
	running_jobs = []
	while x<6.0:
		maintain_queue_size(running_jobs, 4)
		atoms = [utils.Atom('Pb', 0,0,0), utils.Atom('Cl', x,0,0)]
		#name = 'PbCl+_x%.2f' % x
		name = 'PbCl+_x%.2f_mp2' % x
		print 'Running', name
		#job = orca.job(name, '! B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} Grid3 FinalGrid5', atoms, queue=None, grad=True, charge_and_multiplicity='1 1')
		#job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='1 1')
		#running_jobs.append(job)
		#labels for cml file
		for a in atoms:
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')
		x+=0.25

def pairwise_pbcl_qz(): #compare energies for name and old_name below: QZ is apparently not neccessary
	x = 2.0
	running_jobs = []
	while x<6.0:
		maintain_queue_size(running_jobs, 4)
		atoms = [utils.Atom('Pb', 0,0,0), utils.Atom('Cl', x,0,0)]
		name = 'PbCl+_x%.2f_mp2_qz' % x
		old_name = 'PbCl+_x%.2f_mp2' % x
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ Def2-QZVPP(-g,-f) ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='1 1', previous=old_name)
		running_jobs.append(job)
		#labels for cml file
		for a in atoms:
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')
		x+=0.25

def steps_pbcl2():
	running_jobs = []
	x = 0.0
	images = []
	while x<4.0:
		y = 2.0
		while y<6.0:
			maintain_queue_size(running_jobs, 4)
			atoms = [utils.Atom('Pb', 0,0,0), utils.Atom('Cl', 2.5,0,0), utils.Atom('Cl', -x,y,0)]
			name = 'PbCl2_x%.2f_y%.2f_mp2' % (x,y)
			print 'Running', name
			job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True)
			running_jobs.append(job)
			images.append(atoms)
			for a in atoms: #labels for cml file
				if a.element=='Pb': a.label='907'
				if a.element=='Cl': a.label='344'
			files.write_cml(atoms, name='orca/'+name+'/system.cml')
			y+=0.5
		x+=0.5
	files.write_xyz(images)

def analyze_md_states_Pb2Cl3():
	states = files.read_xyz('md_states')[:50]
	running_jobs = []
	for i,atoms in enumerate(states):
		maintain_queue_size(running_jobs, 4)
		name = 'Pb2Cl3+_MD_%d_mp2' % i
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='1 1', previous=('Pb2Cl3+_MD_0_mp2' if i>4 else None))
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

analyze_md_states_Pb2Cl3()


