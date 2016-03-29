from merlin import *
import shutil, re

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
	Pb: utils.Struct(index=Pb, index2=Pb_, element_name='Pb', element=82, mass=207.2, charge=0.4, vdw_e=10.1, vdw_r=3.0),
	I: utils.Struct(index=I, index2=I_, element_name='I', element=53, mass=126.9, charge=-0.2, vdw_e=10.1, vdw_r=2.5),
	Cl: utils.Struct(index=I, index2=Cl_, element_name='Cl', element=17, mass=35.435, charge=-0.2, vdw_e=10.1, vdw_r=2.0),
}

system = utils.System(box_size=[30, 30, 30], name=run_name)
DMSO = utils.Molecule('molecules/dmso')
PbCl2 = utils.Molecule('molecules/PbCl2', extra_parameters=extra, check_charges=False)
system.add(PbCl2)
system.add(PbCl2, 8.0)
#files.packmol(system, (DMSO,PbCl2), (20,1), 0.4, random_seed)
#files.write_xyz(system.atoms)

os.chdir('lammps')
files.write_lammps_data(system)

shutil.copy('md_input.tersoff', system.name+'_input.tersoff')

commands = ('''units real
atom_style full
pair_style hybrid/overlay lj/cut/coul/inout 0.2 0.0 12 tersoff
bond_style harmonic
angle_style harmonic
dihedral_style opls
special_bonds lj/coul 0.0 0.0 0.5

boundary p p p
read_data	'''+system.name+'''.data
''').splitlines()

for line in open(system.name+'_input.tersoff'):
	if line.startswith('# Charges:'): charges = line.split()[2:] #assumes there are 2 tersoff types
	if line.startswith('# LJ-sigma:'): lj_sigma = line.split()[2:]
	if line.startswith('# LJ-epsilon:'): lj_epsilon = line.split()[2:]
	
tersoff_cutoffs = re.findall('\n     +\S+ +\S+ +\S+ +\S+ +(\S+)', open(system.name+'_input.tersoff').read())
inner_cutoffs = [float(tersoff_cutoffs[0]), float(tersoff_cutoffs[-1])] #assumes there are 2 tersoff types

for i in range(len(system.atom_types)):
	for j in range(i, len(system.atom_types)):
		if i>=len(charges) or j>=len(charges):
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f 0.0' % (i+1, j+1, (system.atom_types[i].vdw_e*system.atom_types[j].vdw_e)**0.5, (system.atom_types[i].vdw_r*system.atom_types[j].vdw_r)**0.5) )
		else:
			commands.append('pair_coeff %d %d lj/cut/coul/inout %f %f %f' % (i+1, j+1, (float(lj_epsilon[i])*float(lj_epsilon[j]))**0.5, (float(lj_sigma[i])*float(lj_sigma[j]))**0.5, (inner_cutoffs[i]*inner_cutoffs[j])**0.5) )
			commands.append('set type %d charge %f' % (i+1, float(charges[i])) ) #todo: set inner cutoff dynamically from Tersoff file

lmp = utils.Struct()
lmp.file = open(system.name+'.in', 'w')
def writeline(line):
	lmp.file.write(line+'\n')
lmp.command = writeline
for line in commands:
	lmp.command(line)

#run LAMMPS
lmp.command('pair_coeff * * tersoff '+system.name+'_input.tersoff Pb Cl '+(' NULL'*(len(system.atom_types)-2)) )

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
dump	1 all xyz 100 '''+system.name+'''.xyz
thermo_style custom step temp epair emol vol
thermo 1000
fix motion all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
velocity all create 300.0 '''+random_seed+''' rot yes dist gaussian
timestep 1.0
run 10000
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

