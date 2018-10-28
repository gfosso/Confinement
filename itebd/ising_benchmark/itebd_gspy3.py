""" iTEBD code to find the ground state of 
the 1D Ising model on an infinite chain.
The results are compared to the exact results.
Frank Pollmann, frankp@pks.mpg.de"""

import numpy as np
from scipy import integrate
from scipy.linalg import expm 
from ham import *
from mps import *
import time
start=time.time()
#Hamiltonian

hlong=0.0
htras=0.5
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
	#B.append(np.zeros([2,1,1])); B[-1][0,0,0]=1
	#s.append(np.ones([1]))
        B.append(np.random.rand(2,10,10))
        s.append(np.random.rand(10))

# Generate the two-site time evolution operator
#H_bond = np.array([[J,gx*0.5,gx*0.5,0], [gx*0.5,-J,0,gx*0.5], [gx*0.5,0,-J,gx*0.5], [0,gx*0.5,gx*0.5,J]] )
#H_bond += np.diag([-gz,0,0,gz])
H_bond=Ham
U = np.reshape(expm(-delta*H_bond),(2,2,2,2))

# Perform the imaginary time evolution alternating on A and B bonds
for step in range(0, N):
    s,B = evol(s,B,U,chi,d)

    
# compute magnetization
mag=magnetization(s,B,d)
print("sigmazeta =", mag)


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
end=time.time()
print(end-start)
