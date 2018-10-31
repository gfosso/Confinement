""" iTEBD code to quench the chain"""

import matplotlib.pyplot as plt
from ham import *
from mps import *
import time
start=time.time()
#Hamiltonian
htras=0.25
gort=0.100
gpar=3.0
# diagonal part
Ham = np.diag([ -gort*0.5*SzSz(conf,0,1) -gort*0.5*SzSz(conf,2,3) for conf in range(hilbertsize)])
Ham += np.diag([-gpar*SzSz(conf,0,2) -gpar*SzSz(conf,1,3) for conf in range(hilbertsize)])
# off-diagonal part
for conf in range(hilbertsize):
        value, newconf = Spinflip(conf,0,2)
        Ham[newconf,conf] -=value     
        value, newconf = Spinflip(conf,1,3)
        Ham[newconf,conf] -=value     
#transverse external field
for conf in range(hilbertsize):
        value, newconf = Sx(conf,0)
        Ham[newconf,conf] -= htras*value     
        value, newconf = Sx(conf,1)
        Ham[newconf,conf] -= htras*value     
        value, newconf = Sx(conf,2)
        Ham[newconf,conf] -= htras*value     
        value, newconf = Sx(conf,3)
        Ham[newconf,conf] -= htras*value     

print(Ham)

# First define the parameters of the model / simulation
J=-1.; chi=150; d=4; delta=0.1; T=5; L=int(T/delta);

# Generate the two-site time evolution operator
H_bond = Ham
U = np.reshape(expm(-complex(0,delta)*H_bond),(4,4,4,4))
corr=[]
# Perform the real time evolution alternating on A and B bonds
for step in range(0, L): 
   v=[corrszsz(i,s,B,d) for i in range(0,5)]
   corr.append(v)
 #   corr.append(magnetization(s,B,d))
   s,B=evol(s,B,U,chi,d)


                # compute magnetization
#print "sigmazeta =", np.mean(mag)

t=np.arange(0,T,delta)
#distance= np.arange(0,40,1)
#plt.contourf(distance,t,corr,10,cmap='RdGy')
#plt.colorbar()
#plt.show()

#to save data to a file use
np.savetxt('corr.out',corr)


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
