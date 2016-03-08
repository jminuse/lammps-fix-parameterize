from merlin import *
import shutil

I_ = 66
H_ = 54
N_ = 53
Pb_ = 111

Pb = 907
I = 838

extra = {
	(H_, I_): (100.0, 2.1), 
	(N_, H_, I_): (10.0, 180.0), 
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=0.1, vdw_r=2.0, D0=5.0, alpha=1.5, r0=2.8),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=0.1, vdw_r=2.0, D0=5.0, alpha=1.5, r0=2.8),
	(Pb_, I_): (100.0, 2.9), 
	(I_, Pb_, I_): (10.0, 95.0),
	(13, 53, 54, 66): (0.0,0.0,0.0),
	(54, 53, 54, 66): (0.0,0.0,0.0),
}

system = utils.System(box_size=[1e3, 1e3, 1e3], name='test0')

for root, dirs, file_list in os.walk("data"):
	for ff in file_list:
		if ff.endswith('.cml'):
			total = utils.Molecule('data/'+ff, extra_parameters=extra, check_charges=False)
			system.add(total, len(system.molecules)*200.0)
system.box_size[0] = len(system.molecules)*400+200

for root, dirs, file_list in os.walk('orca'):
	count = 0
	for d in dirs:
		print d
	for ff in file_list:
		if ff.endswith('.out'):
			print ff
				
		#for step in range(20):
		#		name = 'PbI2_r%d' % step
		#		if not name.startswith('PbI') : continue #for PbI testing
		#		if not name.endswith('_def2SVP'): continue
		#		energy, atoms = g09.parse_atoms(name, check_convergence=False)
		#		if len(atoms)>3: continue
		#		if any([utils.dist(atoms[0], a)>3.5 for a in atoms]) and len(atoms)<6: continue
		#		total = utils.Molecule('gaussian/'+name, extra_parameters=extra, check_charges=False)
		#		total.energy = energy*627.509 #convert energy from Hartree to kcal/mol
		#		total.element_string = ' '.join( [a.element for a in total.atoms] )
		#		print total.element_string
		#		for i,a in enumerate(total.atoms):
		#			b = atoms[i]
		#			a.x, a.y, a.z = b.x, b.y, b.z
		#			a.fx, a.fy, a.fz = [f*1185.8113 for f in (b.fx, b.fy, b.fz)] # convert forces from Hartree/Bohr to kcal/mol / Angstrom
		#		system.add(total, count*200.0)
		#		count += 1

exit()

os.chdir('lammps')
files.write_lammps_data(system)

shutil.copy('md_tersoff.tersoff', system.name+'.tersoff')

f = open('target_forces.txt', 'w')
for a in system.atoms:
	a.fx, a.fy, a.fz = 0.0, 0.0, 0.0
	f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
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

compute atom_pe all pe/atom
''').splitlines()
lmp = utils.Struct()
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

#run LAMMPS
seed = sys.argv[1] if len(sys.argv)>1 else '1'
lmp.command('pair_coeff * * tersoff '+system.name+'.tersoff Pb I '+(' NULL'*(len(system.atom_types)-2)) ) #is it possible to do this with the LAMMPS set command?

for t in system.atom_types:
	if hasattr(t,'vdw_e'):
		lmp.command('set type %d charge %f' % (t.lammps_type, t.charge))
		lmp.command('pair_coeff %d * lj/cut/coul/inout %f	%f' % (t.lammps_type, t.vdw_e, t.vdw_r) )
for t in system.bond_types:
	lmp.command('bond_coeff %d	%f %f' % (t.lammps_type, t.e, t.r) )
for t in system.angle_types:
	lmp.command('angle_coeff %d	%f %f' % (t.lammps_type, t.e, t.angle) )
for t in system.dihedral_types:
	lmp.command('dihedral_coeff %d	%f %f %f %f' % ((t.lammps_type,)+t.e))

commands = '''
dump 1 all xyz 10000 '''+system.name+'''.xyz
thermo 100
fix params all parameterize target_forces.txt upper_bounds.tersoff lower_bounds.tersoff '''+seed+'''
run 10000
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('/fs/home/jms875/build/lammps/lammps-7Dec15/src/lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

