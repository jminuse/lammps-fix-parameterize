from merlin import *
import shutil

run_name = sys.argv[1] if len(sys.argv)>1 else 'test'
random_seed = sys.argv[2] if len(sys.argv)>2 else '1'

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
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=0.1, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=0.1, vdw_r=2.5),
	Cl: utils.Struct(index=I, index2=Cl_, element_name='Cl', element=17, mass=35.435, charge=-0.2, vdw_e=0.1, vdw_r=2.0),
}

system = utils.System(box_size=[1e3, 1e3, 1e3], name=run_name)

systems_by_composition = {}

for outer in ['/fs/home/jms875/build/lammps/lammps-7Dec15/src/test/']:
	directories = next(os.walk(outer+'orca'))[1]
	for directory in directories:
		if not os.path.isfile(outer+'orca/'+directory+'/'+directory+'.orca.engrad'): continue
		atoms, energy = orca.engrad_read(outer+'orca/'+directory+'/'+directory+'.orca.engrad')
		if len(atoms)>2: continue
		with_bonds = utils.Molecule(outer+'orca/'+directory+'/system.cml', extra_parameters=extra, check_charges=False)
		for a,b in zip(atoms,with_bonds.atoms):
			convert = 627.51/0.529177249 #Hartee/Bohr to kcal/mole-Angstrom
			b.fx, b.fy, b.fz = a.fx*convert, a.fy*convert, a.fz*convert
		with_bonds.energy = energy
		composition = ' '.join(sorted([a.element for a in atoms]))
		if composition not in systems_by_composition:
			systems_by_composition[composition] = []
		systems_by_composition[composition].append(with_bonds)

for composition in systems_by_composition: #within each type of system, lowest energy must be first and equal to 0.0
	systems_by_composition[composition].sort(key=lambda s:s.energy)
	for s in systems_by_composition[composition]:
		s.energy -= systems_by_composition[composition][0].energy
		system.add(s, len(system.molecules)*1000.0)

#energies.sort()
#print [e-energies[0] for e in energies]
#exit()

system.box_size[0] = len(system.molecules)*400+200

os.chdir('lammps')
files.write_lammps_data(system)
f = open('test_log.txt','w')
for t in system.atom_types:
	f.write(str(t)+' ')
f.close()

shutil.copy('input.tersoff', system.name+'.tersoff')
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
pair_style hybrid/overlay lj/cut/coul/inout 0.0001 3 15 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary p p p
read_data	'''+system.name+'''.data

pair_coeff * * lj/cut/coul/inout 0.0 1.0 0
''').splitlines()
lmp = utils.Struct()
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

#run LAMMPS
lmp.command('pair_coeff * * tersoff '+system.name+'.tersoff Pb Cl '+(' NULL'*(len(system.atom_types)-2)) )

#for t in system.atom_types:
#	if hasattr(t,'vdw_e'):
#		lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
#		lmp.command('pair_coeff %d * lj/cut/coul/inout %f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
for t in system.bond_types:
	lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
for t in system.angle_types:
	lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
for t in system.dihedral_types:
	lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

commands = '''
compute atom_pe all pe/atom
thermo 0
neigh_modify once yes
fix params all parameterize '''+system.name+'''_forces.txt  '''+system.name+'''_energies.txt  '''+system.name+'''_upper.tersoff '''+system.name+'''_lower.tersoff '''+system.name+'''_best.tersoff '''+random_seed+'''
run 100000000
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

