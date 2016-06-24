### AFH 6.22.16
#This file takes the input to run_min_params and creates an upper and lower bounds file for it based on a 50% 
#upper and lower range
#NOTE this only works for Pb Cl systems without H3
###
from merlin import *

def main():

	os.chdir('lammps')

	upper = open('upper_bounds_noH3.tersoff',"r+") #open files for doing MCSMRFF without H3
	myfile = open('input_noH3.tersoff',"r+")
	lower = open('lower_bounds_noH3.tersoff',"r+")

	tmp_upper = open('upper_tmp.txt','w+')
	tmp_lower = open('lower_tmp.txt','w+')
	
	for line in upper:
		if line.startswith('#'):
			tmp_upper.write(line) #write comments to upper

	for line in lower:
		if line.startswith('#'):
			tmp_lower.write(line) #write comments to lower

	for line in myfile: #find upper and lower bounds for params
		if not line.startswith('#'):
			linesplit = line.split()
			if len(linesplit) > 0:
				if (linesplit[0] == 'Pb') or (linesplit[0] == 'Cl'):
					
					atoms = "%s  %s  %s " % (linesplit[0], linesplit[1], linesplit[2])
					m = linesplit[3]
					gamma = linesplit[4]
					lambda3 = linesplit[5]
					c = linesplit[6]
					d = linesplit[7]
					costheta0 = linesplit[8]
				
					m_upper = float(m) #find upper bounds
					gamma_upper = findhigher(gamma)
					lambda3_upper = findhigher(lambda3)
					c_upper = findhigher(c)
					d_upper = findhigher(d)
					costheta0_upper = findhigher(costheta0)

					m_lower = float(m) #find lower bounds
					gamma_lower = findlower(gamma)
					lambda3_lower = findlower(lambda3)
					c_lower = findlower(c)
					d_lower = findlower(d)
					costheta0_lower = findlower(costheta0)
				
					tmp_upper.write("%s	%f   %f   %f   %f   %f   %f\n" % (atoms, m_upper, gamma_upper, lambda3_upper, c_upper, d_upper, costheta0_upper))

					tmp_lower.write("%s	%f   %f   %f   %f   %f   %f\n" % (atoms, m_lower, gamma_lower, lambda3_lower, c_lower, d_lower, costheta0_lower))
				
				else:
					n = linesplit[0] 
					beta = linesplit[1]
					lambda2 = linesplit[2]
					B = linesplit[3]
					R = linesplit[4]
					D = linesplit[5]
					lambda1 = linesplit[6]
					A = linesplit[7]
				 
					n_upper = findhigher(n) #find upper bounds
					beta_upper = findhigher(beta)
					lambda2_upper = findhigher(lambda2)
					B_upper = findhigher(B)
					R_upper = findhigher(R)
					D_upper = findhigher(D)
					lambda1_upper = findhigher(lambda1)
					A_upper = findhigher(A)
					
					n_lower = findlower(n) #find lower bounds
					beta_lower = findlower(beta)
					lambda2_lower = findlower(lambda2)
					B_lower = findlower(B)
					R_lower = findlower(R)
					D_lower = findlower(D)
					lambda1_lower = findlower(lambda1)
					A_lower = findlower(A)

					tmp_upper.write("\t\t%f   %f   %f   %f   %f   %f   %f   %f\n\n" % (n_upper, beta_upper, lambda2_upper, B_upper, R_upper, D_upper, lambda1_upper, A_upper))

					tmp_lower.write("\t\t%f   %f   %f   %f   %f   %f   %f   %f\n\n" % (n_lower, beta_lower, lambda2_lower, B_lower, R_lower, D_lower, lambda1_lower, A_lower))
		
	tmp_upper.close()
	tmp_lower.close()
	myfile.close()
	tmp_upper = open('upper_tmp.txt','r')
	tmp_lower = open('lower_tmp.txt','r')
	upper = open('upper_bounds_noH3.tersoff',"w")
	lower = open('lower_bounds_noH3.tersoff',"w")
	upper.truncate(0)
	lower.truncate(0)
	for line in tmp_upper: #write the data to the correct files
		upper.write(line)
	for line in tmp_lower:
		lower.write(line)

	tmp_upper.close() #close all the files
	tmp_lower.close()
	upper.close()
	lower.close()

def findlower(num):
	num = float(num)
	low = num - num*.5
	return low

def findhigher(num):
	num = float(num)
	high = num + num*.5
	return high	


main()
