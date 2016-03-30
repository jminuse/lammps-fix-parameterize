from merlin import *

import matplotlib.pyplot as plt
import numpy as np
import time, copy, random

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

def analyze_md_states_PbCl2():
	states = files.read_xyz('PbCl2_states')[10::20]
	running_jobs = []
	for i,atoms in enumerate(states):
		maintain_queue_size(running_jobs, 2)
		name = 'PbCl2_MD_%d_mp2' % (i+6)
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='0 1', previous='PbCl2_MD_1_mp2')
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def analyze_md_states_Pb2Cl4():
	states = files.read_xyz('Pb2Cl4_states')[5::5]
	running_jobs = []
	for i,atoms in enumerate(states):
		maintain_queue_size(running_jobs, 2)
		name = 'Pb2Cl4_MD_%d_mp2' % (i)
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue=None, grad=True, charge_and_multiplicity='0 1', previous=('Pb2Cl4_MD_0_mp2' if i>1 else None))
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')
		
def opt_Pb2Cl4():
	#orca.job('Pb2Cl4_MD_opt3_mp2', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv Opt LooseOpt', queue=None, previous='Pb2Cl4_MD_opt_mp2', procs=1, extra_section='%scf SOSCFStart 0.00003 end')
	#orca.job('Pb2Cl4_MD_opt4_mp2', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv Opt LooseOpt', atoms=files.read_xyz('states'), queue=None, previous='Pb2Cl4_MD_opt_mp2', procs=1, extra_section='%scf SOSCFStart 0.00003 end')
	#orca.job('Pb2Cl4_MD_opt6_mp2', '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', queue=None, grad=True, previous='Pb2Cl4_MD_opt4_mp2')
	for name in ['Pb2Cl4_MD_opt5_mp2', 'Pb2Cl4_MD_opt6_mp2', 'Pb2Cl_3+_opt_mp2']:
		atoms = orca.read(name).atoms
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def analyze_md_states_Pb2Cl_3():
	#states = files.read_xyz('Pb2Cl_3+')[5::5]
	states = files.read_xyz('out')[40::40]
	running_jobs = []
	for i,atoms in enumerate(states):
		maintain_queue_size(running_jobs, 4)
		name = 'Pb2Cl_3+_MD_%d_mp2' % (i+33)
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6', atoms, queue=None, grad=True, charge_and_multiplicity='3 1', previous='Pb2Cl_3+_MD_26_mp2')
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def opt_Pb2Cl_3():
	#orca.job('Pb2Cl_3+_opt', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv Opt LooseOpt', queue=None, previous='Pb2Cl_3+_MD_0_mp2', procs=1, extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='3 1')
	orca.job('Pb2Cl_3+_opt_mp2', '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', queue=None, grad=True, previous='Pb2Cl_3+_opt', charge_and_multiplicity='3 1')

def analyze_md_states_Pb4Cl8():
	states = files.read_xyz('out')[:20]
	running_jobs = []
	for i,atoms in enumerate(states):
		maintain_queue_size(running_jobs, 4)
		name = 'Pb4Cl8_MD_%d_mp2' % (i)
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6', atoms, queue=None, grad=True, previous=('Pb4Cl8_MD_0_mp2' if i>3 else None))
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def PbCl6_neg4():
	atoms = files.read_xyz('molecules/PbCl6_4-')
	#orca.job('PbCl6_-4_opt', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv Opt LooseOpt', atoms, queue=None, procs=1, charge_and_multiplicity='-4 1', extra_section='%scf SOSCFStart 0.00003 end')
	orca.job('PbCl6_-4_opt_mp2', '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv GridX6 Opt LooseOpt', queue='batch', grad=True, previous='PbCl6_-4_sp_0_mp2', extra_section='%scf SOSCFStart 0.00033 directresetfreq 1 MaxIter 500 end')

def PbCl6_neg4_states():
	base_atoms = orca.read('PbCl6_-4_opt').atoms
	running_jobs = []
	for i in range(20):
		maintain_queue_size(running_jobs, 4)
		atoms = copy.deepcopy(base_atoms)
		for a in atoms:
			a.x = random.gauss(a.x, i*0.1)
		name = 'PbCl6_-4_sp_%d' % i
		print 'Running', name
		job = orca.job(name, '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv SP', atoms, queue=None, charge_and_multiplicity='-4 1', extra_section='%scf SOSCFStart 0.00003 end', previous='PbCl6_-4_opt')
		running_jobs.append(job)

def show_PbCl6_neg4_states():
	jobs = [orca.read('PbCl6_-4_sp_%d' % i) for i in range(12)]
	jobs.sort(key=lambda j:j.energy)
	print [ 1000*(j.energy-jobs[0].energy) for j in jobs ]
	files.write_xyz([j.atoms for j in jobs])

def PbCl6_neg4_states_mp2():
	running_jobs = []
	for i in range(1,12):
		maintain_queue_size(running_jobs, 4)
		old_name = 'PbCl6_-4_sp_%d' % i
		name = 'PbCl6_-4_sp_%d_mp2' % i
		atoms = orca.read(old_name).atoms
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv GridX6', atoms, queue=None, grad=True, previous='PbCl6_-4_sp_0_mp2', extra_section='%scf SOSCFStart 0.00033 directresetfreq 1 MaxIter 500 end')
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

PbCl6_neg4()
#PbCl6_neg4_states_mp2()

def PbCl6_neg4_states_grad():
	running_jobs = []
	for i in range(12):
		#maintain_queue_size(running_jobs, 4)
		old_name = 'PbCl6_-4_sp_%d' % i
		name = 'PbCl6_-4_sp2_%d' % i
		atoms = orca.read(old_name).atoms
		print 'Running', name
		job = orca.job(name, '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid3 FinalGrid5 SlowConv SP', queue='batch', grad=True, charge_and_multiplicity='-4 1', extra_section='%scf SOSCFStart 0.00003 end', previous=old_name)
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

