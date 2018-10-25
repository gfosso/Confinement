""" iTEBD code to find the ground state of 
the 1D Ising model on an infinite chain.
The results are compared to the exact results.
Frank Pollmann, frankp@pks.mpg.de"""

import numpy as np
from scipy import integrate
from scipy.linalg import expm 
from ham import *

#Hamiltonian

hlong=0.0
htras=0.0
# diagonal part
Ham=[]
Ham = np.diag([ -4.*SzSz(conf,0,1) -hlong*Sz(conf,0) -hlong*Sz(conf,1) for conf in range(hilbertsize)])
# off-diagonal part
for conf in range(hilbertsize):
        value, newconf = Sx(conf,0)
        Ham[newconf,conf] -= htras*value     
        value, newconf = Sx(conf,1)
        Ham[newconf,conf] -= htras*value     
print(Ham)


# First define the parameters of the model / simulation
J=-1.; chi=100; d=2; delta=0.01; N=3000;
B=[];s=[]
for i in range(2):
#	B.append(np.zeros([2,1,1])); B[-1][0,0,0]=1
#	s.append(np.ones([1]))
        B.append(np.random.rand(2,10,10))
        s.append(np.random.rand(10))

# Generate the two-site time evolution operator
#H_bond = np.array([[J,gx*0.5,gx*0.5,0], [gx*0.5,-J,0,gx*0.5], [gx*0.5,0,-J,gx*0.5], [0,gx*0.5,gx*0.5,J]] )
#H_bond += np.diag([-gz,0,0,gz])
H_bond=Ham
U = np.reshape(expm(-delta*H_bond),(2,2,2,2))

# Perform the imaginary time evolution alternating on A and B bonds
for step in range(0, N):
	for i_bond in [0,1]:
		ia = np.mod(i_bond-1,2); ib = np.mod(i_bond,2); ic = np.mod(i_bond+1,2)
		chia = B[ib].shape[1]; chic = B[ic].shape[2]

		# Construct theta matrix and time evolution #
		theta = np.tensordot(B[ib],B[ic],axes=(2,1)) # i a j b
		theta = np.tensordot(U,theta,axes=([2,3],[0,2])) # ip jp a b 
		theta = np.tensordot(np.diag(s[ia]),theta,axes=([1,2])) # a ip jp b 
		theta = np.reshape(np.transpose(theta,(1,0,2,3)),(d*chia,d*chic)) # ip a jp b

		# Schmidt decomposition #
		X, Y, Z = np.linalg.svd(theta,full_matrices=0)
		chi2 = np.min([np.sum(Y>10.**(-10)), chi])	
		piv = np.zeros(len(Y), np.bool)
		piv[(np.argsort(Y)[::-1])[:chi2]] = True

		Y = Y[piv]; invsq = np.sqrt(sum(Y**2))
		X = X[:,piv] 
		Z = Z[piv,:]
		
		# Obtain the new values for B and s #
		s[ib] = Y/invsq 
		X=np.reshape(X,(d,chia,chi2))
		X = np.transpose(np.tensordot(np.diag(s[ia]**(-1)),X,axes=(1,1)),(1,0,2))
		B[ib] = np.tensordot(X, np.diag(s[ib]),axes=(2,0))
		
		B[ic] = np.transpose(np.reshape(Z,(chi2,d,chic)),(1,0,2))

# compute magnetization
sz=np.diag([Sz(conf,0) for conf in range(0,2)])
mag=[]
for i_bond in range(2):
	sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
	C = np.tensordot(sB,sz,axes=(1,0))
	sB=np.conj(sB)
	mag.append(np.squeeze(np.tensordot(sB,C,axes=([0,2,1],[0,1,2]))).item()) 
print("sigmazeta =", np.mean(mag))


# Get the bond energies
E=[]
for i_bond in range(2):
    BB = np.tensordot(B[i_bond],B[np.mod(i_bond+1,2)],axes=(2,1))
    sBB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),BB,axes=(1,1))
    C = np.tensordot(sBB,np.reshape(H_bond,[d,d,d,d]),axes=([1,2],[2,3]))
    sBB=np.conj(sBB)
    E.append(np.squeeze(np.tensordot(sBB,C,axes=([0,3,1,2],[0,1,2,3]))).item()) 
print("E_iTEBD =",np.mean(E))

f = lambda k,g : -2*np.sqrt(1+g**2-2*g*np.cos(k))/np.pi/2.
E0_exact = integrate.quad(f, 0, np.pi, args=(htras,))[0]
print("E_exact =", E0_exact)
