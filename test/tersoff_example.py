from merlin import *

def read_tersoff_file(filename):
	# format of a single entry (one or more lines):
	#   element 1, element 2, element 3, 
	#           m, gamma, lambda3, c, d, costheta0, 
	#           n, beta, lambda2, B, R, D, lambda1, A
	tersoff_params = []
	contents = ''.join([line for line in open(filename) if not line.startswith('#')])
	lines = contents.split('\n\n') #expects entries to be separated by a blank line
	for line in lines:
		line = line.replace('\n','') #an entry can span any number of lines
		columns = line.split()
		if len(columns)==17:
			e1, e2, e3, m, gamma, lambda3, c, d, costheta0, n, beta, lambda2, B, R, D, lambda1, A = [ (float(s) if s[-1].isdigit() else s) for s in columns]
			tersoff_params.append( utils.Struct(e1=e1, e2=e2, e3=e3, m=m, gamma=gamma, lambda3=lambda3, c=c, d=d, costheta0=costheta0, n=n, beta=beta, lambda2=lambda2, B=B, R=R, D=D, lambda1=lambda1, A=A) )
	
	return tersoff_params

def write_tersoff_file(system, best=False,error=-1):
	if best:
		f = open(system.name+'_best.tersoff', 'w')
		f.write('# Error: ' + str(error)+'\n')
	else:
		f = open(system.name+'.tersoff', 'w')

	pb_i_i = [t for t in system.tersoff_params if (t.e1,t.e2,t.e3)==('Pb','I','I')][0]
	for i,t in enumerate(system.tersoff_params):
		if (t.e1,t.e2,t.e3)==('I','Pb','Pb'):
			n = t
			t = pb_i_i
			f.write(('%-2s '*3+('%12.6g '*6)+'\n         '+('%12.6g '*8)+'\n\n') % (n.e1, n.e2, n.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
		else:
			f.write(('%-2s '*3+('%12.6g '*6)+'\n         '+('%12.6g '*8)+'\n\n') % (t.e1, t.e2, t.e3, t.m, t.gamma, t.lambda3, t.c, t.d, t.costheta0, t.n, t.beta, t.lambda2, t.B, t.R, t.D, t.lambda1, t.A))
	f.close()

def set_lammps_parameters(system):
	write_tersoff_file(system)
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

def run(system):
	#run LAMMPS
	seed = str(1)
	set_lammps_parameters(system)
	
	f = open('target_forces.txt', 'w')
	for a in system.atoms:
		a.fx, a.fy, a.fz = 0.0, 0.0, 0.0
		f.write("%e\n%e\n%e\n" % (a.fx, a.fy, a.fz) )
	f.close()
	
	commands = '''
dump 1 all xyz 10000 '''+system.name+'''.xyz
thermo 100
fix params all parameterize target_forces.txt '''+seed+'''
run 10000
	'''
	for line in commands.splitlines():
		lmp.command(line)
	
	lmp.file.close()
	os.system('/fs/home/jms875/build/lammps/lammps-7Dec15/src/lmp_serial -in %s.in -log %s.log' % (system.name,system.name))

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

system = utils.System(box_size=[20, 20, 20], name='md_tersoff')

m = utils.Molecule('PbI2', extra_parameters=extra, check_charges=False)
for x in range(-1,1,1):
	for y in range(0,1,1):
		for z in range(0,1,1):
			system.add(m, x*6, y*6, z*4)

os.chdir('lammps')
files.write_lammps_data(system)

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

system.tersoff_params = read_tersoff_file('input.tersoff')

run(system)

