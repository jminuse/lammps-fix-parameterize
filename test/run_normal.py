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
acetone = utils.Molecule('molecules/acetone')
PbCl2 = utils.Molecule('molecules/PbCl2', extra_parameters=extra, check_charges=False)
MACl = utils.Molecule('molecules/MACl', extra_parameters=extra, check_charges=False)

if 0:
	PbCl6_neg4 = utils.Molecule('molecules/PbCl6_4-', extra_parameters=extra, check_charges=False)
	system.add(PbCl6_neg4)
elif 0: #solvent+solute cluster
	for xi in range(-1,1):
		for yi in range(-1,1):
			for zi in range(1):
				system.add(PbCl2, xi*6, yi*6, zi*6)
elif 0: #Pb2Cl
	Pb2Cl_3 = utils.Molecule('molecules/Pb2Cl_3+', extra_parameters=extra, check_charges=False)
	system.add(Pb2Cl_3)
elif 0: #solvent+solute cluster
	for xi in range(3):
		for yi in range(3):
			system.add(DMSO, (xi-0.5)*6, (yi-0.5)*6)
	pb_ion = utils.Molecule('molecules/pb2+', extra_parameters=extra, check_charges=False)
	system.add(pb_ion, 0, 0, 5)
	system.add(PbCl2, 10, 10, 10)
elif 0: #MACl+PbCl2
	for xi in range(3):
		for yi in range(3):
			system.add(MACl, (xi-0.5)*10, (yi-0.5)*10)
	#system.add(PbCl2, 0, 0, 5)
elif 1: #perovskite cubic structure
	PbMACl3 = utils.Molecule('molecules/PbMACl3', extra_parameters=extra, check_charges=False)
	L = 6.6
	N = 1
	system.box_size=[N*L, N*L, N*L]
	for xi in range(N):
		for yi in range(N):
			for zi in range(N):
				x, y, z = (xi-0.5)*L, (yi-0.5)*L, (zi-0.5)*L
				system.add(PbMACl3, x, y, z)
else: #packed, solvated system
	files.packmol(system, (DMSO,acetone,PbCl2), (5,5,1), 1.5, random_seed)

files.write_xyz(system.atoms) #for testing purposes
#exit()

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
neigh_modify every 1 check yes delay 0
dump	1 all xyz 100 '''+system.name+'''.xyz
thermo_style custom step temp epair emol vol
thermo 1000
minimize 0.0 1.0e-8 1000 100000
min_style fire
minimize 0.0 1.0e-8 1000 100000
#fix motion all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
fix motion all nvt temp 400.0 400.0 100.0
velocity all create 400.0 '''+random_seed+''' rot yes dist gaussian
timestep 1.0
run 10000
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

