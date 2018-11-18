import numpy as np
#maronna
#size
L=10
#hilbertsize
hilbertsize=2**L
delta=-13322

def binconf(c): return np.binary_repr(c,L)

def readsite(conf,i): return (conf&(1<<i))>>i

def SzSz(conf,i,j):
    si = readsite(conf,i)-1/2
    sj = readsite(conf,j)-1/2
    return si*sj

def Sz(conf,i):
    return readsite(conf,i)-0.5

def XXZHam(conf):
    return sum([- delta*( 4.*SzSz(conf,i,(i+1)%L) + 1. ) for i in range(L) ])

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


def count(conf):
	return sum(readsite(conf,i) for i in range(L))
#for conf in range(hilbertsize):
#        print(binconf(conf),binconf(Spinflip(conf,0,1)[1]))


#reduced hilbert space in tot_sz=0 naive way

c=[]
for conf in range(hilbertsize):
	if count(conf) == 0.5*L:
		c.append(conf)




# Hamiltonian
# diagonal part
Ham = np.diag([XXZHam(c[i]) for i in range(len(c))])
# off-diagonal part
for j in range(len(c)):
    for i in range(L):
        value, newconf = Spinflip(c[j],i,(i+1)%L)
        Ham[c.index(newconf),j] -= 2.*value     
        
        
print(Ham)
en ,pin = np.linalg.eigh(Ham)
print(en[0]/L)



