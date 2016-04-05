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
	orca.job('PbCl6_-4_opt_mp2', '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv Opt LooseOpt', queue=None, extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1', previous='PbCl6_-4_opt')

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
		#maintain_queue_size(running_jobs, 4)
		old_name = 'PbCl6_-4_sp_%d' % i
		name = 'PbCl6_-4_sp_%d_mp2' % i
		atoms = orca.read(old_name).atoms
		print 'Running', name
		job = orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue='batch', grad=True, previous='PbCl6_-4_opt_mp2', extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1')
		running_jobs.append(job)
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

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

def PbCl6_neg4_md_mp2():
	states = files.read_xyz('pbcl6_states')[1::10]
	for i,atoms in enumerate(states):
		name = 'PbCl6_-4_md_%d_mp2' % i
		orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue='batch', grad=True, previous='PbCl6_-4_opt_mp2', extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1')
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

#PbCl6_neg4_md_mp2()
#PbCl6_neg4_states_mp2()
#PbCl6_neg4()

def PbCl6_neg4_opt_sp_mp2():
	states = orca.read('PbCl6_-4_opt_mp2').frames
	for i,atoms in enumerate(states):
		name = 'PbCl6_-4_opt_%d_mp2' % i
		orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', atoms, queue='batch', grad=True, previous='PbCl6_-4_opt_mp2', extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1')
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def PbCl6_bigger_basis():
	directories = next(os.walk('orca'))[1]
	for directory in directories:
		name = directory
		if 'PbCl6' not in name: continue
		if not os.path.isfile('orca/'+name+'/'+name+'.orca.engrad'): continue
		if not os.path.isfile('orca/'+name+'/system.cml'): continue
		try:
			atoms, energy = orca.engrad_read('orca/'+name+'/'+name+'.orca.engrad', pos='Ang')
		except IndexError:
			continue
		if name.endswith('_ma'): continue
		print name
		old_name = name
		name = old_name + '_ma' # Add minimally-augmented basis set, because anion. Without this, anion errors can be quite large. 
		#orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ ma-def2-TZVP ECP{def2-TZVP} AutoAux TIGHTSCF Grid5 FinalGrid6 SlowConv', queue='batch', grad=True, previous=old_name, extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1')
		atoms = orca.read(old_name).atoms
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def PbCl6_bigger_bigger_basis():
	directories = next(os.walk('orca'))[1]
	for directory in directories:
		name = directory
		if 'PbCl6' not in name or 'mp2' not in name: continue
		if not os.path.isfile('orca/'+name+'/'+name+'.orca.engrad'): continue
		if not os.path.isfile('orca/'+name+'/system.cml'): continue
		try:
			atoms, energy = orca.engrad_read('orca/'+name+'/'+name+'.orca.engrad', pos='Ang')
		except IndexError:
			continue
		if not name.endswith('_ma'): continue
		#print name
		old_name = name
		name = old_name + '2' # From TZ to QZ
		orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ ma-def2-QZVPP ECP{def2-SD} AutoAux TIGHTSCF Grid5 FinalGrid6 SlowConv', queue='batch', grad=True, previous=old_name, extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1')
		atoms = orca.read(old_name).atoms
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def PbCl6_bigger_bigger_basis_more_memory():
	directories = next(os.walk('orca'))[1]
	for directory in directories:
		name = directory
		if not name.endswith('_ma2'): continue
		if not os.path.exists('orca/'+name+'/'+name+'.orca.gbw'): continue
		#print name
		old_name = name
		name = old_name[:-1] + '3' # more memory
		if os.path.exists('orca/'+name): continue
		orca.job(name, '! RIJCOSX RI-B2PLYP D3BJ ma-def2-QZVPP ECP{def2-SD} AutoAux TIGHTSCF Grid5 FinalGrid6 SlowConv', queue='batch', grad=True, previous=old_name, extra_section='%scf SOSCFStart 0.00003 end', charge_and_multiplicity='-4 1', mem=500)
		atoms = orca.read(old_name[:-1]).atoms
		for a in atoms: #labels for cml file
			if a.element=='Pb': a.label='907'
			if a.element=='Cl': a.label='344'
		files.write_cml(atoms, name='orca/'+name+'/system.cml')

def compare_TZ_QZ_PbCl6_anion():
	directories = next(os.walk('orca'))[1]
	jobs = []
	for directory in directories:
		name = directory
		if name.endswith('_ma'):
			try:
				ma = orca.read(name)
				ma3 = orca.read(name+'3')
				jobs.append( (ma3.energy,ma.energy) )
			except:
				print name, 'failed'
	jobs.sort()
	for j in jobs:
		e_qz, e_tz = j[0]-jobs[0][0], j[1]-jobs[0][1]
		print e_qz*1000, e_tz*1000
	
def opt_unit_cell():
	def initial():
		atoms = utils.Molecule('molecules/PbMACl3', check_charges=False).atoms
		orca.job('PbMACl3_opt', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} Grid3 FinalGrid5 Opt LooseOpt', atoms, queue=None)
	orca.job('PbMACl3_opt2', '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} Grid3 FinalGrid5 Opt LooseOpt', queue=None, previous='PbMACl3_opt', extra_section='%geom Constraints '+(' '.join(['{C %d C}'%i for i in range(1,9)]))+' end end')

def mod_unit_cell():
	atoms = orca.read('PbMACl3_opt').atoms
	for i in range(10,15):
		for a in atoms[8:]:
			a.x = random.gauss(a.x, (i+1)*0.1)
		orca.job('PbMACl3_sp_%d'%i, '! B97-D3 GCP(DFT/TZ) def2-TZVP ECP{def2-TZVP} Grid3 FinalGrid5 SP', atoms, queue=None, previous='PbMACl3_opt').wait()

def mp2_unit_cell():
	for i in range(15):
		#orca.job('PbMACl3_mp2_%d'%i, '! RIJCOSX RI-B2PLYP D3BJ def2-TZVP ECP{def2-TZVP} TIGHTSCF Grid5 FinalGrid6 SlowConv', queue='batch', grad=True, previous='PbMACl3_sp_%d'%i)
		system = utils.Molecule('molecules/PbMACl3', check_charges=False)
		atoms = orca.read('PbMACl3_sp_%d'%i).atoms
		for a,b in zip(system.atoms,atoms):
			a.x,a.y,a.z = b.x,b.y,b.z
		files.write_cml(system, 'orca/PbMACl3_mp2_%d/system.cml'%i)

def blaire_dimer_jobs():
	d = '/fs/home/bas348/perovskites/xyz/'
	names = ['2pbcl3ma_p', '2pbcl3ma', '2pbcl3ma_4dmso_p', '2pbcl3ma_3dmso_p', '2pbcl3ma_2dmso_p', '2pbcl3ma_1dmso_p', '2pbcl3ma_4dmso', '2pbcl3ma_3dmso', '2pbcl3ma_3dmso', '2pbcl3ma_2dmso', '2pbcl3ma_1dmso']
	for name in names:
		name = '2pbcl3ma_1dmso'
		atoms = files.read_xyz(d+name+'.xyz')
		#orca.job(name+'_opt', '! Opt B97-D3 def2-TZVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO', atoms, queue='batch', extra_section='%cosmo  SMD true  solvent "DMSO"  end')
		
		
		extra_basis = '%basis\n'
		elements = dict([(a.element,True) for a in atoms]).keys()
		for e in elements:
			if e=='Pb':
				extra_basis += 'NewECP '+e+' "def2-SD" "def2-TZVP" end NewAuxGTO '+e+' "def2-TZVP/J" end\n'
			elif e=='S':
				extra_basis += 'newGTO '+e+' "def2-TZVP" end NewAuxGTO '+e+' "def2-TZVP/J" end\n'
			else:
				extra_basis += 'newGTO '+e+' "def2-SVP" end NewAuxGTO '+e+' "def2-SVP/J" end\n'
		extra_basis += 'end'
		

		orca.job(name+'_opt_dz', '! Opt B97-D3 def2-SVP ECP{def2-TZVP} COSMO printbasis', atoms, queue=None, extra_section='%cosmo  SMD true  solvent "DMSO"  end '+extra_basis).wait()
		exit()

def read_ssh():
	os.system("ssh vesuvius 'cd /tmp/icse_1465927.1; tail *.out;'")

def new_dimer_jobs(): #apart energy = -6105.575540556, joined energy = -6105.570657715
	names = ['2pbcl3ma_apart_5dmso', '2pbcl3ma_joined_5dmso', '2pbcl3ma_joined_5dmso_2']
	for name in names:
		atoms = files.read_xyz('molecules/'+name+'.xyz')
		
		for i in range(3):
			for a in atoms:
				for b in atoms:
					if a is not b:
						if utils.dist(a,b)<1e-4:
							atoms.remove(a)
							continue
		
		extra_basis = '%basis\n'
		elements = dict([(a.element,True) for a in atoms]).keys()
		for e in elements:
			if e=='Pb':
				extra_basis += 'NewECP '+e+' "def2-SD" "def2-TZVP" end NewAuxGTO '+e+' "def2-TZVP/J" end\n'
			elif e=='S':
				extra_basis += 'newGTO '+e+' "def2-TZVP" end NewAuxGTO '+e+' "def2-TZVP/J" end\n'
			else:
				extra_basis += 'newGTO '+e+' "def2-SVP" end NewAuxGTO '+e+' "def2-SVP/J" end\n'
		extra_basis += 'end'
		
		orca.job(name+'_opt_dz', '! Opt B97-D3 def2-SVP ECP{def2-TZVP} COSMO printbasis', atoms, queue='batch', extra_section='%cosmo  SMD true  solvent "DMSO"  end\n'+extra_basis)
	
def test_orca_basis_syntax():
	names = ['2pbcl3ma_apart_5dmso', '2pbcl3ma_joined_5dmso', '2pbcl3ma_joined_5dmso_2']
	for name in names:
		atoms = files.read_xyz('molecules/'+name+'.xyz')
		orca.job(name+'_opt', '! Opt B97-D3 def2-SVP GCP(DFT/TZ) ECP{def2-TZVP} COSMO printbasis', atoms, queue='batch', extra_section='%cosmo  SMD true  solvent "DMSO"  end\n%basis aux auto NewECP Pb "def2-SD" "def2-TZVP" end NewGTO S "def2-TZVP" end end')



