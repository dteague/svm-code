#!/bin/bash
for i in 0 1 2 3
do
    tmp=$[$i**2 + 10]
    echo $[2*$i*$tmp/3 + $tmp]
done
