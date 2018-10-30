import matplotlib.pyplot as plt
import numpy as np

T=30
dist=30
delta=0.1
t=np.arange(0,T,delta)
distance= np.arange(0,dist,1)
plt.contourf(distance,t,np.loadtxt('corr.out'),20)
        #cmap='RdGy')
plt.colorbar()
plt.show()
