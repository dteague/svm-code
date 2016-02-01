import os
import subprocess
import time

i = 13
x = 1
change = 5
i_max = i+change #18
j_max = 4
front = 0

wavel = 100
width = 30



filetext = "width		:" + str(width) + "\n"
filetext += "lambda		:" + str(wavel) + "d-9\n"
filetext += "#_timesteps	:3000\n"
filetext += "dt		:0.05\n"
filetext += "Z_charge	:2\n"
filetext += "T_of_output	:500\n"
filetext += "half_lat_num	:100\n"
filetext += "lat_width	:.02\n"
filetext += "taylor		:.false.\n"
filetext += "disp_wavefunc	:.true.\n"
filetext += "disp_density	:.true.\n"
filetext += "is_helium	:.true.\n"
filetext += "is_hydrogen	=.false.\n"
filetext += "I_0\t:"

while i < i_max:
	j = 0
	while j < j_max:
		if j == 0:
			front = 10
		elif j==1:
			front = 18
		elif j==2:
			front = 30
		elif j==3:
			front = 56
		input = open('input.inp','w')
		input.write(filetext)
		input.write(str(front)+"d"+str(i))
		input.close()
		os.system('../../thecode')
		gs_file = str(front) + 'd' + str(i) + '.dat'
		os.system('mv groundstate.dat ' + gs_file)
		
		os.system('tail -n 250 ' + gs_file + ' > ' + str(x) + '.dat')
		
		x = x + 1
		j = j + 1
	i=i+1
max=j_max*change
file_out = 'knee_f' + str(wavel)+ '_w' + str(width)
	
	
proc = subprocess.Popen(['gnuplot','-p'], 
                        shell=True,
                        stdin=subprocess.PIPE,
                        )


proc.stdin.write('f(x)=a\n')
proc.stdin.write('set print "' + file_out + '.dat"\n')
i=1
while i < (1+max):
	proc.stdin.write('fit f(x) "%d.dat" via a\n' % i)
	proc.stdin.write('print %d,a\n' % i)
	i +=1
proc.stdin.write('set print\n')	
proc.stdin.write('quit\n')

time.sleep(60)
os.system('mkdir ' + file_out)
os.system('mv ' + file_out + '.dat ./' + file_out)
os.system('mv nohup.out ./' + file_out)
os.system('mv fit.log ./' + file_out)
i=0
while i < 10:
	os.system('mv ' + str(i) + '* ./' + file_out)
	i+=1
