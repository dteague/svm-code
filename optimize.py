import os


x0_start=0.5
x0_step=0.5
a0_start=2
a0_step=0.1

temp=""
i=0
space_limit = 1
j=0
status = open('status.out','w')
while j < space_limit:
	x0 = x0_start+x0_step*j
	status.write("\t")
	status.write(str(x0))
	j=j+1

status.write("\n")
while i < 1:
	j=0
	a0 = a0_start+a0_step*i
	status.write(str(a0))
	while j < space_limit:
		x0 = x0_start+x0_step*j
		
		fileText = "4\n"
		fileText += "0\n"
		fileText += "0\n"	
		fileText += str(x0) + "\n"
		fileText += str(a0) + "\n"
		
		basisInp = open('basis.inp', 'w')
		basisInp.write(fileText)
		basisInp.close()
	
	
		
		os.system('./basis_o_w')
		os.system('./thecode')
		
		file = open('eigenenergies.dat', 'r')
		energy = file.readline()
		
		status.write("\t")
		if energy == temp:
			status.write("NaN")
		else:
			e_str=str(energy)
			status.write(e_str[4:10])
		temp = energy
		j=j+1
	i=i+1
	status.write("\n")
