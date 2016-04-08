from merlin import *
import shutil, hashlib, re

run_name = sys.argv[1] if len(sys.argv)>1 else 'test'
random_seed = sys.argv[2] if len(sys.argv)>2 else str(int(hashlib.md5(run_name).hexdigest(), 16)%(2**16))

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
	#system.xlo,system.xhi = -N*L*0.5, N*L*0.5
	#system.ylo,system.yhi = -N*L*0.5, N*L*0.5
	#system.zlo,system.zhi = -N*L*0.5, N*L*0.5
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
timestep 0.5
run 10000
'''
for line in commands.splitlines():
	lmp.command(line)

lmp.file.close()
os.system('../../lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

