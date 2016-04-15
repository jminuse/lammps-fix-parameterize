from merlin import *
import shutil, hashlib, re

run_name = sys.argv[1] if len(sys.argv)>1 else 'test'
optimization_method = sys.argv[2] if len(sys.argv)>2 else 'SBPLX'
random_seed = int(sys.argv[3]) if len(sys.argv)>3 else int(hashlib.md5(run_name).hexdigest(), 16)%(2**16)

I_ = 66
Cl_ = 21
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838
Cl = 344
HN = 233

extra = {
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=10.1, vdw_r=3.0),
}

system = utils.System(box_size=[1e3, 100.0, 100.0], name=run_name)

systems_by_composition = {}

for outer in ['/fs/home/jms875/build/lammps/lammps-7Dec15/src/test/']:
	directories = next(os.walk(outer+'orca'))[1]
	for directory in directories:
		name = directory
		if not os.path.isfile(outer+'orca/'+name+'/'+name+'.orca.engrad'): continue
		if not os.path.isfile(outer+'orca/'+name+'/system.cml'): continue
		try:
			atoms, energy = orca.engrad_read(outer+'orca/'+name+'/'+name+'.orca.engrad', pos='Ang')
		except IndexError:
			continue
		if all([s not in name for s in ['PbMACl3_mp2_', 'PbCl6NH3CH3-2_mp2_sp_']]): continue  #'PbMACl3x2_mp2_', 'PbMACl3_mp2_', 
		#	if len(atoms)!=3 or 'mp2' not in name or 'qz' in name or 'PbCl2_' not in name: continue
		if '_md_' in name: continue
			#if len(atoms)<4 or len(atoms)>6 or len(atoms)==5 or 'mp2' not in name or 'qz' in name: continue
		#if '-4' in name and not name.endswith('ma3'): continue # strong anion without augmented basis
		#if 'PbCl6_' in name and not ('_ma3' in name and '_opt_' in name): continue
		with_bonds = utils.Molecule(outer+'orca/'+name+'/system.cml', extra_parameters=extra, check_charges=False)
		for b in with_bonds.bonds:
			r = utils.dist(*b.atoms)
			if r > 2.5:
				raise Exception('Bond too large (%.5g A) in %s' % (r,name))
		for a,b in zip(atoms,with_bonds.atoms):
			convert = 627.51/0.529177249 #Hartee/Bohr to kcal/mole-Angstrom
			b.fx, b.fy, b.fz = a.fx*convert, a.fy*convert, a.fz*convert
			if utils.dist(a,b)>1e-4:
				raise Exception('Atoms too different:', (a.x,a.y,a.z), (b.x,b.y,b.z))
		with_bonds.energy = energy
		with_bonds.name = name
		composition = ' '.join(sorted([a.element for a in atoms]))
		if composition not in systems_by_composition:
			systems_by_composition[composition] = []
		systems_by_composition[composition].append(with_bonds)

xyz_atoms = []

for composition in systems_by_composition: #within each type of system, lowest energy must be first and equal to 0.0
	systems_by_composition[composition].sort(key=lambda s:s.energy)
	baseline_energy = systems_by_composition[composition][0].energy
	for s in systems_by_composition[composition]:
		s.energy -= baseline_energy
		s.energy *= 627.5 #Convert Hartree to kcal/mol
		if s.energy > 500.0: continue #don't use high-energy systems, because these will not likely be sampled in MD
		print composition, s.name, s.energy #for testing purposes
		xyz_atoms.append(s.atoms) #for testing purposes
		system.add(s, len(system.molecules)*1000.0)

system.xhi = len(system.molecules)*1000.0+100.0

files.write_xyz(xyz_atoms, 'states')
#exit()

os.chdir('lammps')
files.write_lammps_data(system)
f = open(system.name+'_data.txt', 'w')
for composition in systems_by_composition:
	for s in systems_by_composition[composition]:
		f.write(s.name+'\n')
f.close()

shutil.copy('input.tersoff', system.name+'_input.tersoff')
shutil.copy('upper_bounds.tersoff', system.name+'_upper.tersoff')
shutil.copy('lower_bounds.tersoff', system.name+'_lower.tersoff')

# write forces to a file
f = open(system.name+'_forces.txt', 'w')
for a in system.atoms:
	f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
f.close()
# write energies to a file
f = open(system.name+'_energies.txt', 'w')
for m in system.molecules:
	f.write("%e\n" % (m.energy) )
f.close()

commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/inout 0.2 3.5 15 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary f f f
read_data	'''+system.name+'''.data
''').splitlines()

tersoff_types = [t for t in system.atom_types if t.index in [Pb,Cl,HN]]

for line in open(system.name+'_input.tersoff'):
	if line.startswith('# Charges:'): charges = [float(x) for x in line.split()[2:]]
	if line.startswith('# LJ-sigma:'): lj_sigma = [float(x) for x in line.split()[2:]]
	if line.startswith('# LJ-epsilon:'): lj_epsilon = [float(x) for x in line.split()[2:]]

tersoff_strings = re.findall('\n'+('(\S+) +'*9)[:-2]+' *\n +'+('(\S+) +'*8)[:-2], open(system.name+'_input.tersoff').read())
inner_cutoffs = {}
for type_i in tersoff_types:
	for type_j in tersoff_types:
		for s in tersoff_strings:
			types = s[:3]
			R, D = float(s[13]), float(s[14])
			if types == (type_i.element_name, type_j.element_name, type_j.element_name):
				inner_cutoffs[ (type_i,type_j) ] = R+D

#for indices,cutoff in inner_cutoffs.items():
#	print indices[0].element_name, indices[1].element_name, cutoff

index = 0
for t in system.atom_types:
	if t in tersoff_types:
		t.charge = charges[index]
		t.vdw_e = lj_epsilon[index]
		t.vdw_r = lj_sigma[index]
		index += 1

for i in range(len(system.atom_types)):
	for j in range(i, len(system.atom_types)):
		type_i = system.atom_types[i]
		type_j = system.atom_types[j]
		commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (i+1, j+1, (type_i.vdw_e*type_j.vdw_e)**0.5, (type_i.vdw_r*type_j.vdw_r)**0.5, inner_cutoffs[ (i,j) ] if (i,j) in inner_cutoffs else 0.0) )
	commands.append('set type %d charge %f' % (i+1, type_i.charge) )

lmp = utils.Struct()
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

lmp.command('pair_coeff * * tersoff '+system.name+'_input.tersoff '+(' '.join([ (t.element_name if t in tersoff_types else 'NULL') for t in system.atom_types])))

for t in system.bond_types:
	lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
for t in system.angle_types:
	lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
for t in system.dihedral_types:
	lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

commands = '''
compute atom_pe all pe/atom
compute sum_pe all reduce sum c_atom_pe
neigh_modify once yes

min_style params
min_modify '''+run_name+''' '''+optimization_method+''' '''+str(random_seed)+'''
minimize 0.01 0.01 1728000000 1728000000 #with 224 atoms, does 2e4 steps/second. One day = 1728000000, 40% of 2^32. 
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

