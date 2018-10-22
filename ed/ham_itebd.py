import numpy as np

#size
L=2
#hilbertsize
hilbertsize=2**(2*L)
gort=0.5
gpar=0.5

def binconf(c): return np.binary_repr(c,L)

def readsite(conf,i): return (conf&(1<<i))>>i

def SzSz(conf,i,j):
    si = readsite(conf,i)-1/2
    sj = readsite(conf,j)-1/2
    return si*sj

def Sz(conf,i):
    return readsite(conf,i)-0.5


#for conf in range(hilbertsize):
#    print(IsingHam(conf))

def Sx(conf,i):
    return 0.5, conf^(1<<i)


#for conf in range(hilbertsize):
#    print(binconf(conf),binconf(Sx(conf,0)[1]))


def Spinflip(conf,i,j):
    "SxSx + SySy = 0.5(S+S- + S-S+) term"
    if readsite(conf,i) != readsite(conf,j):
        return 0.5, conf^((1<<i)^(1<<j))
    else: return 0.0, conf


#for conf in range(hilbertsize):
#        print(binconf(conf),binconf(Spinflip(conf,0,1)[1]))




#Hamiltonian
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


