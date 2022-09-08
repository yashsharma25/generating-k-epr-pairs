# example of "4-pairable" graph state:
# one can use LOCC to create any pattern of two EPR states on any subset of 4 qubits

import numpy as np
import itertools

####################### helper functions #############################

def null2(A):
# Argument: binary matrix A
# Returns a matrix whose rows span the nullspace of A over the binary field GF(2)  
    rows,n = A.shape
    X = np.identity(n,dtype=int)
    for i in range(rows):
        y = np.dot(A[i,:], X) % 2
        not_y = (y + 1) % 2
        good = X[:,np.nonzero(not_y)]
        good = good[:,0,:]
        bad = X[:, np.nonzero(y)]
        bad = bad[:,0,:]
        if bad.shape[1]>0 :
            bad = np.add(bad,  np.roll(bad, 1, axis=1) )  % 2
            bad = np.delete(bad, 0, axis=1)
            X = np.concatenate((good, bad), axis=1)
    return np.transpose(X)

def rank2(A):
	KerA = null2(A)
	return A.shape[1]-KerA.shape[0]

def test_epr(G,a,b):
# Argument: stabilizer tableaux of size s x 2n
# each row is a stabilizer, columns range(n) = X part, columns range(n,2n) = Z part
# checks whether qubits a and b are maximally entangled with each other
	Sa = rank2(G[:,[a,a+n]])
	Sb = rank2(G[:,[b,b+n]])
	Sab = rank2(G[:,[a,b,a+n,b+n]])
	if Sa==2 and Sb==2 and Sab==2:
		return True
	else:
		return False

####################################################################

# number of qubits
n = 10

# define adjacency matrix of the graph state
A = np.zeros((n,n),dtype=int)
for i in range(n):
	j = (i+1) % n
	A[i,j] = 1
	A[j,i] = 1
	j = (i+int(n/2)) % n
	A[i,j] = 1
	A[j,i] = 1

# compute all pairings of the set {0,1,...,n-1}
S4 = itertools.permutations(range(4))
four_elem_pairings = [sigma for sigma in S4 if sigma[0]<sigma[1] and sigma[2]<sigma[3] and sigma[0]<sigma[2] ]
n_elem_pairings = []
for c in itertools.combinations(range(n), 4):
	for p in four_elem_pairings:
		n_elem_pairings.append([(c[p[0]],c[p[1]]),(c[p[2]],c[p[3]])])
assert(len(n_elem_pairings)==int(3*n*(n-1)*(n-2)*(n-3)/24))


# compute a list of all Pauli measurement bases
xyz = []
for q in range(n-4):
	xyz.append(['Z','X','Y'])
pauli_basis = []
for elem in itertools.product(*xyz):
	pauli_basis.append(elem)


# stabilizer tableaux of the graph state
G = np.zeros((n,2*n),dtype=int)
for i in range(n):
	G[i,i] = 1
	G[i,n:2*n] = A[i,:] 
GT = np.transpose(G)

# symplectric inner product matrix
Lambda = np.zeros((2*n,2*n),dtype=int)
for q in range(n):
	Lambda[q][q+n] = 1
	Lambda[q+n][q] = 1

for pp in n_elem_pairings:
	isOK = False
	for meas in pauli_basis:
		a,b = pp[0]
		c,d = pp[1]
		qubits_to_measure = np.sort(list(set(range(n))-set([a,b,c,d])))

		# stabilizer tableaux for the measured qubits
		GM = np.zeros((n-4,2*n),dtype=int)
		for i in range(n-4):
			q = qubits_to_measure[i]
			if meas[i]=='X':
				GM[i,q] = 1
			if meas[i]=='Y':
				GM[i,q] = 1
				GM[i,q+n] = 1
			if meas[i]=='Z':
				GM[i,q+n] = 1

		# find graph state stabilizers commuting with all measurements
		Gcom = np.dot(null2(np.dot(GM,np.dot(Lambda,GT)) % 2),G) % 2

		if test_epr(Gcom,a,b) and test_epr(Gcom,c,d):
			print('EPR pairs',pp,'measurement basis=',end='')
			cnt = 0
			for q in range(n):
				if q in qubits_to_measure:
					print(meas[cnt],sep='',end='')
					cnt = cnt + 1
				else:
					print('*',sep='',end='')
			print('')
			isOK = True
			break
	if not(isOK):
		print('EPR pairs ',pp,' cannot be generated')
		print('Terminating the search...') 
		exit()
