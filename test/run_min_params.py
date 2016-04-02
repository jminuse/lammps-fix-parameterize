from merlin import *
import shutil, hashlib, re

run_name = sys.argv[1] if len(sys.argv)>1 else 'test'
optimization_method = sys.argv[2] if len(sys.argv)>2 else 'SBPLX'
random_seed = int(hashlib.md5(run_name).hexdigest(), 16)%(2**16)

I_ = 66
Cl_ = 21
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838
Cl = 344

extra = {
	(H_, I_): (100.0, 2.1), 
	(N_, H_, I_): (10.0, 180.0), 
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=10.1, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=10.1, vdw_r=2.5),
	Cl: utils.Struct(index=I, index2=Cl_, element_name='Cl', element=17, mass=35.435, charge=-0.2, vdw_e=10.1, vdw_r=2.0),
}

system = utils.System(box_size=[1e3, 50, 50], name=run_name)

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
		#if 'PbMACl3_mp2_' not in name: continue
		if len(atoms)>6 or 'mp2' not in name or 'qz' in name or len(atoms)==5: continue
		#if '-4' in name and not name.endswith('ma3'): continue # strong anion without augmented basis
		#if 'PbCl6_' in name and not ('_ma3' in name and '_opt_' in name): continue
		with_bonds = utils.Molecule(outer+'orca/'+name+'/system.cml', extra_parameters=extra, check_charges=False)
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

system.box_size[0] = len(system.molecules)*1000.0*2+200.0
count = 0
for m in system.molecules:
	files.write_xyz(m.atoms,str(count))
	count += 1

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


for line in open(system.name+'_input.tersoff'):
	if line.startswith('# Charges:'): charges = line.split()[2:] #assumes there are 2 tersoff types
	if line.startswith('# LJ-sigma:'): lj_sigma = line.split()[2:]
	if line.startswith('# LJ-epsilon:'): lj_epsilon = line.split()[2:]
	
tersoff_cutoffs = re.findall('\n     +\S+ +\S+ +\S+ +\S+ +(\S+)', open(system.name+'_input.tersoff').read())
inner_cutoffs = [float(tersoff_cutoffs[0]), float(tersoff_cutoffs[-1])] #assumes there are 2 tersoff types

pb_type = 0
cl_type = 1

for i in range(len(system.atom_types)):
	for j in range(i, len(system.atom_types)):
		if i>=len(charges) or j>=len(charges):
			if i==pb_type:
				commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f 0.0' % (i+1, j+1, (system.atom_types[i].vdw_e*system.atom_types[j].vdw_e)**0.5, (system.atom_types[i].vdw_r*system.atom_types[j].vdw_r)**0.5-0.0) )
			elif i==cl_type:
				commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f 0.0' % (i+1, j+1, (system.atom_types[i].vdw_e*system.atom_types[j].vdw_e)**0.5, (system.atom_types[i].vdw_r*system.atom_types[j].vdw_r)**0.5+0.0) )
			else:
				commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f 0.0' % (i+1, j+1, (system.atom_types[i].vdw_e*system.atom_types[j].vdw_e)**0.5, (system.atom_types[i].vdw_r*system.atom_types[j].vdw_r)**0.5) )
		else:
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (i+1, j+1, (float(lj_epsilon[i])*float(lj_epsilon[j]))**0.5, (float(lj_sigma[i])*float(lj_sigma[j]))**0.5, (inner_cutoffs[i]*inner_cutoffs[j])**0.5) )
			commands.append('set type %d charge %f' % (i+1, float(charges[i])) )

lmp = utils.Struct()
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

lmp.command('pair_coeff * * tersoff '+system.name+'_input.tersoff Pb Cl '+(' NULL'*(len(system.atom_types)-2)) )

for t in system.bond_types:
	lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
for t in system.angle_types:
	lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
for t in system.dihedral_types:
	lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

commands = '''
compute atom_pe all pe/atom
compute sum_pe all reduce sum c_atom_pe
#thermo_style custom c_sum_pe
#thermo 1
neigh_modify once yes

min_style params
min_modify '''+run_name+''' '''+optimization_method+''' '''+str(random_seed)+'''
minimize 0.01 0.01 1728000000 1728000000 #with 224 atoms, does 2e4 steps/second. One day = 1728000000, 40% of 2^32. 
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

