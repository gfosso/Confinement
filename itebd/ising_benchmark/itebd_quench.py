""" iTEBD code to quench the chain"""
import matplotlib.pyplot as plt
from ham import *
from correlators import *
import time
start=time.time()
#Hamiltonian

hlong=0.2
htras=0.25
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
J=-1.; chi=100; d=2; delta=0.1; T=20; L=int(T/delta);
#B1[0]=B[0].astype(complex);l1[0]=l[0].astype(complex)

sz=np.diag([Sz(conf,0) for conf in range(0,2)])
#sz=np.array([[0,1],[1,0]])

# Generate the two-site time evolution operator
H_bond = Ham
U = np.reshape(expm(-complex(0,delta)*H_bond),(2,2,2,2))
corr=[]
# Perform the real time evolution alternating on A and B bonds
for step in range(0, L): 
    v=[corrszsz(i,s,B) for i in range(1,40)]
    corr.append(v)
#    corr.append(magnetization(s,B))
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
#print "sigmazeta =", np.mean(mag)

time=np.arange(0,T,delta)
distance= np.arange(1,40,1)
plt.contourf(distance,time,corr,20,cmap='RdGy')
plt.colorbar()
plt.show()
# Get the bond energies
#E=[]
#for i_bond in range(2):
#    BB = np.tensordot(B[i_bond],B[np.mod(i_bond+1,2)],axes=(2,1))
#    sBB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),BB,axes=(1,1))
#    C = np.tensordot(sBB,np.reshape(H_bond,[d,d,d,d]),axes=([1,2],[2,3]))
#    sBB=np.conj(sBB)
#    E.append(np.squeeze(np.tensordot(sBB,C,axes=([0,3,1,2],[0,1,2,3]))).item()) 
#print "E_iTEBD =", np.mean(E)

#f = lambda k,g : -2*np.sqrt(1+g**2-2*g*np.cos(k))/np.pi/2.
#E0_exact = integrate.quad(f, 0, np.pi, args=(g,))[0]
#print "E_exact =", E0_exact
end=time.time()
print(end-start)
