### AFH 6.22.16
#This file takes the curent md_input.tersoff file and generates a random set of outputs for a given variable in the
#tersoff functional form.  It then uses run_normal to run these inputs in lammps. You can then analyse the outputs 
#of those runs to determine if any have parameters that are useable
###
from merlin import *
import random
def main():
	i = 0
	while i < 100:
		os.chdir('lammps')

		tmp = open('temp.txt',"w+") #open a temp file 
		myfile = open('md_input.tersoff',"r+")
		
		for line in myfile:
			if not line.startswith('#'):
				linesplit = line.split()
				if len(linesplit) > 0:
					if (linesplit[0] == 'Pb') or (linesplit[0] == 'Cl'):
						
						atoms = "%s  %s  %s " % (linesplit[0], linesplit[1], linesplit[2])
						m = linesplit[3] #get params
						gamma = linesplit[4]
						lambda3 = linesplit[5]
						c = linesplit[6]
						d = linesplit[7]
						costheta0 = linesplit[8]
					
						m_new = float(m) #find random params for input
						gamma_new = random.uniform(findlower(gamma),findhigher(gamma))
						lambda3_new = random.uniform(findlower(lambda3),findhigher(lambda3))
						c_new = random.uniform(findlower(c),findhigher(c))
						d_new = random.uniform(findlower(d),findhigher(d))
						costheta0_new = random.uniform(findlower(costheta0),findhigher(costheta0))
					
						tmp.write("%s	%f   %f   %f   %f   %f   %f\n" % (atoms, m_new, gamma_new, lambda3_new, c_new, d_new, costheta0_new))
					
					else:
						n = linesplit[0] #get params
						beta = linesplit[1]
						lambda2 = linesplit[2]
						B = linesplit[3]
						R = linesplit[4]
						D = linesplit[5]
						lambda1 = linesplit[6]
						A = linesplit[7]
					
						n_new = random.uniform(findlower(n),findhigher(n)) #find random params for input
						beta_new = random.uniform(findlower(beta),findhigher(beta))
						lambda2_new = random.uniform(findlower(lambda2),findhigher(lambda2))
						B_new = random.uniform(findlower(B),findhigher(B))
						R_new = random.uniform(findlower(R),findhigher(R))
						D_new = random.uniform(findlower(D),findhigher(D))
						lambda1_new = random.uniform(findlower(lambda1),findhigher(lambda1))
						A_new = random.uniform(findlower(A),findhigher(A))

						tmp.write("\t\t%f   %f   %f   %f   %f   %f   %f   %f\n\n" % (n_new, beta_new, lambda2_new, B_new, R_new, D_new, lambda1_new, A_new))
			else:
				tmp.write(line)
		tmp.close()
		myfile.close()
		tmp = open('temp.txt','r')
		myfile = open('md_input.tersoff','w')
		myfile.truncate(0)
		for line in tmp: #write the file to md_input.tersoff
			myfile.write(line)
		tmp.close()
		myfile.close()
		os.chdir('..')
		try:
			os.system('python run_normal.py test_%d' % i)
		except:
			pass
		i += 1
	
	

def findlower(num): #find the lower bound
	num = float(num)
	low = num - num*.5
	return low

def findhigher(num): #find the upper bound
	num = float(num)
	high = num + num*.5
	return high


main()
