import numpy as np

#library where the function suitable for building the hamiltonian
#are defined


#size
L=2
#hilbertsize
hilbertsize=2**(L)

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


#IsingHam
def IsingHam(hlong,htras):
    # diagonal part
    Ham=[]
    Ham = np.diag([ -4.*SzSz(conf,0,1) -hlong*Sz(conf,0) -hlong*Sz(conf,1) for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
        value, newconf = Sx(conf,0)
        Ham[newconf,conf] -= htras*value     
        value, newconf = Sx(conf,1)
        Ham[newconf,conf] -= htras*value    
        
    return Ham



