basis = open('basis.dat','r')
output = open('svm.out','w')

nGauss = 0
nGauss =int(basis.readline())
dim = nGauss(nGauss+1)/2

write.output(dim + "\n")
nu = []
i=0
while i < nGauss:
	nu.append(basis.readline())
	lmValues = basis.readline()
	lmValues = ' '.join(lmValues.split())
	position = basis.readline()
	position = ' '.join(position.split())
	for i in position:
		if i != 0:
			exit(1)
print nu	
  
