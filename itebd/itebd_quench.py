""" iTEBD code to quench the chain"""

from ham import *

#Hamiltonian

gort=-0.1
gpar=2.0
# diagonal part
Ham = np.diag([ -gort*0.5*SzSz(conf,0,1) -gort*0.5*SzSz(conf,2,3) for conf in range(hilbertsize)])
Ham += np.diag([-gpar*SzSz(conf,0,2) -gpar*SzSz(conf,1,3) for conf in range(hilbertsize)])
# off-diagonal part
for conf in range(hilbertsize):
        value, newconf = Spinflip(conf,0,2)
        Ham[newconf,conf] -=value     
        value, newconf = Spinflip(conf,1,3)
        Ham[newconf,conf] -=value     
print(Ham)

# First define the parameters of the model / simulation
J=-1.; chi=100; d=4; delta=0.01; T=20; L=int(T/delta);
mag1=[]
#B1[0]=B[0].astype(complex);l1[0]=l[0].astype(complex)

sz=np.diag([Sz(conf,0) for conf in range(0,4)])
# Generate the two-site time evolution operator
H_bond = Ham
U = np.reshape(expm(-complex(0,delta)*H_bond),(4,4,4,4))

# Perform the real time evolution alternating on A and B bonds
for step in range(0, L): 
    mag=[]
    for i_bond in range(2):
        sB = np.tensordot(np.diag(s[np.mod(i_bond-1,2)]),B[i_bond],axes=(1,1))
        C = np.tensordot(sB,sz,axes=(1,0))
        sB=np.conj(sB)
        mag.append(np.squeeze(np.tensordot(sB,C,axes=([0,2,1],[0,1,2]))).item())  
    mag1.append(np.mean(mag))
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