import numpy as np
#maronna
#size
L=4
#hilbertsize
hilbertsize=2**L
delta=-12

def binconf(c): return np.binary_repr(c,L)

def readsite(conf,i): return (conf&(1<<i))>>i

def SzSz(conf,i,j):
    si = readsite(conf,i)-1/2
    sj = readsite(conf,j)-1/2
    return si*sj

def Sz(conf,i):
    return readsite(conf,i)-0.5

def XXZHam(conf):
    return sum([- delta*( 4.*SzSz(conf,i,(i+1)%L)  ) for i in range(L) ])

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

def translate(conf):
	zero=readsite(conf,0)
	return (conf>>1)|(zero<<(L-1))

#gives the lowest integer representative and the periodicity
def lowestrepr(conf):
	conf0=conf
	conf1=conf
	for i in range(L):
		conf=translate(conf)
		if conf1>conf: conf1=conf
		elif conf0==conf: return conf1,i+1 
	return conf1,L

#gives the lowest integer representative and how many translations you need to reach that configuration	
def repr(conf):
	if lowestrepr(conf)[0]==conf: return conf,0
	lowest=lowestrepr(conf)[0]
	for i in range(L):
		conf=translate(conf)
		if conf==lowest:return lowest,i+1

#return True if the state doesn't exist, False if it's already in the configuration vector
def checkstate(lw):
	for i in range(len(c)):
		if c[i]==lw: return False
	return True

#reduced hilbert space in tot_sz=0 naive way
c=[]
for conf in range(hilbertsize):
	if (count(conf) == 0.5*L)&checkstate(lowestrepr(conf)):
		lw=lowestrepr(conf)
		c.append(lw)




# Hamiltonian
# diagonal part
Ham = np.diag([XXZHam(c[i][0]) for i in range(len(c))])
# off-diagonal part
for j in range(len(c)):
    for i in range(L):
        value, newconf = Spinflip(c[j],i,(i+1)%L)
        Ham[c.index(lowestrepr(newconf)),j] -= 2.*value     
        
        
print(Ham)
en ,pin = np.linalg.eigh(Ham)
print(en[0]/L)



