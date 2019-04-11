# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:25:19 2017

@author: James
"""

import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from qutip import *
qutip.settings.has_mkl=False

#Units of frequnecies are in MZ
#Time units will be in us
Deltag = -2*sc.pi*200 #detuning of 493 laser
Deltar = +2*sc.pi*50 #detuning of 650 laser
Om12 = 2*sc.pi*51.5 #Rabi frequrency of |1> to |2> 
Om23 = 2*sc.pi*35.2 #Rabi frequrency of |2> to |3> 
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

# initialise plot and line
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'b-') 

while Deltag < 2*sc.pi*200:
    #Hamiltonian of system RWA
    H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt = expect(sig22, final_state) #Gives expectation value of solved Hamiltonian for excited state
    yg.append(fexpt)
    xg.append(Deltag/(2*sc.pi))
    Deltag += 2*sc.pi*2

plt.plot(xg, yg,'c-')
plt.xlabel('Δg [MHz]')
plt.ylabel('Population in |2>')
plt.show()

Deltar =  -2*sc.pi*200
while Deltar < 2*sc.pi*200:
    #Hamiltonian of system RWA
    Deltag =  -2*sc.pi*20
    H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt = expect(sig22, final_state)  #Gives expectation value of solved Hamiltonian for excited state
    #print(final_state)
    yr.append(fexpt)
    xr.append(Deltar/(2*sc.pi))
    Deltar += 2*sc.pi*2
    
plt.plot(xr, yr, 'r-')
plt.xlabel('Δr [MHz]')
plt.ylabel('Population in |2>')
plt.show()