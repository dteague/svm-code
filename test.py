import subprocess
import sys
proc = subprocess.Popen(['gnuplot','-p'], 
                        shell=True,
                        stdin=subprocess.PIPE,
                        )


proc.stdin.write('f(x)=a\n')
proc.stdin.write('set print "knee75.dat"\n')
i=1
while i < 3:
	proc.stdin.write('fit f(x) "%d.dat" via a\n' % i)
	proc.stdin.write('print %d,a\n' % i)
	i +=1
