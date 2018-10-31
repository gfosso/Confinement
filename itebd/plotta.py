import matplotlib.pyplot as plt
import numpy as np

T=10
dist=10
delta=0.1
t=np.arange(0,T,delta)
distance= np.arange(0,dist,1)
plt.contourf(distance,t,np.loadtxt('corr.out'),20)
        #cmap='RdGy')
plt.title("correlator with gort=0.1  gpar=3 and htras=0.5 to 0.25")
plt.colorbar()
plt.show()
