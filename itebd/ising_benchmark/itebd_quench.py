#!/usr/bin/python3
""" iTEBD code to quench the chain"""

import sys, getopt
import numpy as np
import matplotlib.pyplot as plt
from ham import *
from mps import *
from itebd_gspy3 import *
import time
from scipy.fftpack import fft
start=time.time()

def main(argv):
   # hlong=''
   # outputfile=''
   # try:
   #     opts,args = getopt.getopt(argv,"hg:o:",["gort=","outputfile="])
   # except getopt.GetoptError:
   #     print("itebd_gspy3.py -g <gort value> -o <outputfile>")
   #     sys.exit(2)
   # for opt, arg in opts:
   #     if opt == "-h":
   #         print("itebd_gspy3.py -g <gort value> -o <outputfile>")
   #         sys.exit()
   #     elif opt in ("-g","--gort"):
   #         hlong = float(arg)
   #     elif opt in ("-o","--outputfile"):
   #         outputfile = arg
    hlong=0.
    s,B = ground_state(hlong) 

    #Hamiltonian
    htras=0.25
    hlong=0.2

    # First define the parameters of the model / simulation
    J=-1.; chi=300; d=2; delta=0.01; T=30; L=int(T/delta);

    # Generate the two-site time evolution operator
    H_bond = IsingHam(hlong,htras)
    U = np.reshape(expm(-complex(0,delta)*H_bond),(d,d,d,d))
    corr=[]
    # Perform the real time evolution alternating on A and B bonds
    for step in range(0, L): 
#       v=[corrszsz(i,s,B,d) for i in range(0,10)]
#       corr.append(v)
       corr.append(magnetization(s,B,d))
       s,B=evol(s,B,U,chi,d)
    
    specc=fft(corr)
    plt.plot(np.array([2*np.pi*i/T for i in range(1,L//2)]),np.abs(specc[1:L//2]))
    plt.yscale("log")
    plt.show()
                    # compute magnetization
    #print "sigmazeta =", np.mean(mag)

    #t=np.arange(0,T,delta)
    #distance= np.arange(0,40,1)
    #plt.contourf(distance,t,corr,10,cmap='RdGy')
    #plt.colorbar()
    #plt.show()

    #to save data to a file use
    #np.savetxt(outputfile,corr)


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

if __name__ == "__main__":
    main(sys.argv[1:])
