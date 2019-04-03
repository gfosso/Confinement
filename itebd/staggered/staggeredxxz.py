from mps import *
import sys
from scipy.linalg import expm
import matplotlib.pyplot as plt
from scipy.fftpack import fft

def main(argv):
    cristo=mps(bond_dimension=400)

    #ground state
    hilbertsize=2**2
#    epsilon=0.001
#    gpar=-epsilon**(-1)
#    gpar=-100
#    h=0.01
#
#    HamA = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)-1.*h*Sz(conf,0)+1.*h*Sz(conf,1)   for conf in range(hilbertsize)])
#    HamB = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)+1.*h*Sz(conf,0)-1.*h*Sz(conf,1)   for conf in range(hilbertsize)])
#    # off-diagonal part
#    for conf in range(hilbertsize):
#            value, newconf = Spinflip(conf,0,1)
#            HamA[newconf,conf] -=2.*value     
#            HamB[newconf,conf] -=2.*value     
#    HamA=epsilon*HamA
#    HamB=epsilon*HamB
#    delta=0.001; N=1000000
#    cristo.random_state()
#    UA=np.reshape(expm(-delta*HamA),(2,2,2,2))
#    UB=np.reshape(expm(-delta*HamB),(2,2,2,2))
#    for step in range(N):
#       cristo.evol(UA,UB)

    cristo.product_state(neel=True)

    #quench
    gpar=-5;epsilon=np.abs(1/gpar);h=0.00
    HamA = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)-1.*h*Sz(conf,0)+1.*h*Sz(conf,1)   for conf in range(hilbertsize)])
    HamB = np.diag([ -0.5*gpar*(4.*SzSz(conf,0,1)+1)+1.*h*Sz(conf,0)-1.*h*Sz(conf,1)   for conf in range(hilbertsize)])
    # off-diagonal part
    for conf in range(hilbertsize):
            value, newconf = Spinflip(conf,0,1)
            HamA[newconf,conf] -=2.*value     
            HamB[newconf,conf] -=2.*value     
#    HamA=epsilon*HamA
#    HamB=epsilon*HamB

    delta=0.01;T=7;L=int(T/delta)
    UA=np.reshape(expm(-complex(0,delta)*HamA),(2,2,2,2))
    UB=np.reshape(expm(-complex(0,delta)*HamB),(2,2,2,2))
#    print(U)
    corr=[]
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    corr.append(np.zeros(L))
    spectrum=[]
    spectrum.append(np.zeros(L))
    for step in range(L):
#        corr[0][step]=4.*cristo.expectation_SzSz(1,connected=False)
#        corr[1][step]=4.*cristo.expectation_SzSz(2,connected=False)
#        corr[2][step]=4.*cristo.expectation_SzSz(3,connected=False) 
        corr[3][step]=2.*cristo.staggered_Sz()
        spectrum[0][step]=-sum(cristo.spectrum()*np.log(cristo.spectrum()))
        cristo.evol(UA,UB)
#    np.savetxt("cristo.out",corr)
#    corrf= fft(corr[2])
    time=np.arange(0,T,delta)
#    plt.plot(time,corr[0])
#    plt.plot(time,corr[1])
#    plt.plot(time,corr[2])
#    plt.plot(time,corr[3])
    plt.plot(time,spectrum[0],'r')
##u    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[0])[1:L//2]))
#u    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[1])[1:L//2]))
#u    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[2])[1:L//2]))
#u    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(corr[3])[1:L//2]))
 #   plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(fft(spectrum[0])[1:L//2]))
#    plt.yscale("log")
##   plt.xscale("log")
#    plt.vlines(x=[2.08626,2.39109,2.619213,2.823744],ymin=0,ymax=1000,colors='r')
#    plt.vlines(x=[2.05076115, 2.35786424,2.59221242,2.79920963],ymin=0,ymax=1000,colors='g')
    plt.show()
if __name__=="__main__":
    main(sys.argv[1:])
