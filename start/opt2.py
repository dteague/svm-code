import os


n = 3
ll = 1

x0_start=0.1
x0_step=0.09
x0_n = 10

a0_start=1.1
a0_step=0.2
a0_n = 10

temp=""

last_e = 100
min_e = [0,0,0]
error = 0.0001

fileText = str(n)+"\n"+str(ll)+"\n0\n"


status = open('status.out','w')
for j in range(x0_n):
	x0 = x0_start+x0_step*j
	status.write("\t")
	status.write(str(x0))

status.write("\n")

while abs(last_e - min_e[0]) > error:
    last_e = min_e[0]
    for i in range(a0_n):
        e_prev = 0
            
	a0 = a0_start+a0_step*i
	status.write(str(a0))
	for j in range(x0_n):

            x0 = x0_start+x0_step*j
            
            basisInp = open('basis.inp', 'w')
            basisInp.write(fileText)
            basisInp.write(str(x0) + "\n"+str(a0) + "\n")
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
                status.write(e_str[2:8])
            temp = energy
            if energy != "NaN":
                e_flo = float(energy)

            os.system('rm *.dat')
            if e_prev - e_flo < 0:
                break
            else:
                e_prev = e_flo
            if min_e[0] > e_flo:
                min_e = [e_flo,x0,a0]
                
        status.write("\n")
        
    status.write(str(min_e) + "\n\n")    
    a0_start = min_e[2] - a0_step
    a0_step = 2.*a0_step/a0_n
    x0_start = min_e[1] - x0_step
    x0_step = 2.*x0_step/x0_n
print min_e

basisInp = open('basis.inp', 'w')
basisInp.write(fileText)
basisInp.write(str(min_e[1]) + "\n"+str(min_e[2]) + "\n")
basisInp.close()
                
os.system('./basis_o_w')
os.system('./thecode')
