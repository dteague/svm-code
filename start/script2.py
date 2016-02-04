import os
import subprocess
import time

i_b = 13
i_end = 18

x = 1
j_max = 4
front = [10,18,30,56]

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

for i in range(i_b,i_end+1):
	for j in range(j_max):
		input = open('input.inp','w')
		input.write(filetext)
		input.write(str(front[j])+"d"+str(i))
		input.close()
		os.system('./thecode')
		gs_file = str(front[j]) + 'd' + str(i) + '.dat'
		os.system('mv groundstate.dat ' + gs_file)
		
		os.system('tail -n 250 ' + gs_file + ' > ' + str(x) + '.dat')
		
		x = x + 1

file_out = 'knee_f' + str(wavel)+ '_w' + str(width)
out = open(file_out + ".dat", "w")

for i in range(1,x):
    f = open(str(i)+".dat", "r")
    tot = 0
    for line in f:
        tot += float(line)
    tot = tot/250
    out.write(str(tot) + "\n")

os.system('mkdir ' + file_out)
os.system('mv ' + file_out + '.dat ./' + file_out)

i=0
while i < 10:
	os.system('mv ' + str(i) + '* ./' + file_out)
	i+=1
