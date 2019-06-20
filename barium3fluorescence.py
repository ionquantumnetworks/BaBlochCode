# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:25:19 2017

@author: James
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import scipy.constants as sc
import numpy as np
from qutip import *
qutip.settings.has_mkl=False

#Units of frequnecies are in MZ
#Time units will be in us
Deltag = -2*sc.pi*0 #detuning of 493 laser
Deltar = +2*sc.pi*0 #detuning of 650 laser
Om12 = 2*sc.pi*15.1 #Rabi frequrency of |1> to |2> 
Om23 = 2*sc.pi*5.3 #Rabi frequrency of |2> to |3> 
gamma21 =  2*sc.pi*15.1 #Decay rate of |2> to |1>
gamma23 =  2*sc.pi*5.3 #Decay rate of |2> to |3>

#Basis states for three levels
one, two, three = qutrit_basis()

#Operators between |n> and |m> 
sig11 = one * one.dag()
sig22 = two * two.dag()
sig33 = three * three.dag()
sig12 = one * two.dag()
sig32 = three * two.dag()
sig21 = two * one.dag()
sig23 = two * three.dag()
sig13 = one * three.dag()
e_ops = [sig11,sig22,sig33] #operators to input into mesolve

#Initial state of system
psi0 = one
#Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
c_ops = [C21,C23]

yr = []
xr = []
yg = []
xg = []
x = []
y = []

#2D plots:
# initialise plot and line
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'b-') 

Deltag =  -2*sc.pi*50
while Deltag < 2*sc.pi*50:
    Deltar = 2*sc.pi*8.5
    #Hamiltonian of system RWA
    H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt = expect(sig22, final_state) #Gives expectation value of solved Hamiltonian for excited state
    yg.append(fexpt)
    xg.append(Deltag/(2*sc.pi))
    Deltag += 2*sc.pi*0.5

plt.plot(xg, yg,'c-')
plt.xlabel('Δg [MHz]')
plt.ylabel('Population in |P-state>')
plt.show()

Deltar =  -2*sc.pi*50
while Deltar < 2*sc.pi*50:
    #Hamiltonian of system RWA
    Deltag =  -2*sc.pi*0
    H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt = expect(sig22, final_state)  #Gives expectation value of solved Hamiltonian for excited state
    #print(final_state)
    yr.append(fexpt)
    xr.append(Deltar/(2*sc.pi))
    Deltar += 2*sc.pi*0.5
    
plt.plot(xr, yr, 'r-')
plt.xlabel('Δr [MHz]')
plt.ylabel('Population in |P-state>')
plt.show()

"""
#Attempt to plot 3D plots
fig = plt.figure()
ax = fig.gca(projection='3d')
Z = []
DR = []
DG = []
Deltag = -2*sc.pi*30
while Deltag < 2*sc.pi*30:
    Deltar = -2*sc.pi*50
    while Deltar < 2*sc.pi*50:
        H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
        final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
        fexpt = expect(sig22, final_state)  #Gives expectation value of solved Hamiltonian for excited state
        Z.append(fexpt)
        DR.append(Deltar/(2*sc.pi))
        DG.append(Deltag/(2*sc.pi))
        Deltar += 2*sc.pi*10
    Deltag += 2*sc.pi*10

DR, DG = np.meshgrid(DR, DG)
print(len(DG))
print(len(DR))
print(len(Z))

# Plot the surface.
surf = ax.plot_surface(DR, DG, Z, cmap=cm.cool,linewidth=0, antialiased=True)

# Customize the z axis.
ax.set_zlim(-0.01, 0.15)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
print(Z)
Z = np.array(Z)
zcont = Z.reshape(len(DR))
print(zcont)
plt.figure()
cp = plt.contourf(DR, DG, zcont)
plt.colorbar(cp)
plt.title('Filled Contours Plot')
plt.xlabel('x (cm)')
plt.ylabel('y (cm)')
plt.show()
"""