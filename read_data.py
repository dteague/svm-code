input = open('30d15.dat','r')
output = open('33.dat','w')
i=1
j=1
while i < 1750:
	input.readline()
	i += 1
for line in input:
        output.write(line)
input.close()
output.close()
