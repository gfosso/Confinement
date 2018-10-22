from quantumising import *

L, h = 8, 10

HamH = GeneralOp([ [ -h, Sx(i)] for i in range(L) ])
HamIsing = GeneralOp([ [ -1.0, SzSz(i,(i+1)%L)] for i in range(L) ])
Ham = HamIsing + HamH
print(Ham)

hilbert = IsingHilbert(L=L)
psi0 = Ham.Groundstate(hilbert)

print(psi0.measure(Ham), psi0.measure(Sx(0)), psi0.measure(SzSz(0,L//2)))

