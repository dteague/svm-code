#!/bin/bash

if [ ! -f basis.inp ]
then
    echo "basis.inp doesn't exist, run opt2.py"
    exit
fi

if [ ! -f basis.dat ]
then
   ./basis_o_w
fi

width="30"
lambda="100"

start=18
stop=18
intensity="10 18 30 56"

data="data_wl${width}_w${width}"
mkdir $data

filetext="width\t\t: $width\n"
filetext+="lambda\t\t:${lambda}d-9\n"
filetext+="#_timesteps\t:3000\n"
filetext+="dt\t\t:0.05\n"
filetext+="Z_charge\t:2\n"
filetext+="T_of_output\t:500\n"
filetext+="half_lat_num\t:100\n"
filetext+="lat_width\t:.02\n"
filetext+="taylor\t\t:.false.\n"
filetext+="disp_wavefunc\t:.true.\n"
filetext+="disp_density\t:.true.\n"
filetext+="is_helium\t:.true.\n"
filetext+="is_hydrogen\t=.false."

i_add="I_0\t\t:"

echo -e $filetext > input.inp

i=$start
x=1

while [ $i -le $stop ]
do

    for j in $intensity
    do
	dir="${j}d${i}"
	mkdir $dir
	cat input.inp > ./$dir/input.inp
	echo -e $i_add+$dir >> ./$dir/input.inp
	cp basis.inp $dir
	cp basis.dat $dir
	cd $dir
	../thecode
	tail -n 250 groundstate.dat > ../${data}/${x}.dat
	mv groundstate.dat ../${data}/${j}d${i}.dat
	x=$[$x+1]
	cd ..
	rm -r $dir
    done
    i=$[$i+1]
done
