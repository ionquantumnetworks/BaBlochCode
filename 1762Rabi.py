# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:36:59 2020

@author: James
"""

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
import scipy.constants as sc

#Units
MHz = 1*10**6
kHz = 1*10**3
ns = 1*10**-9
us = 1*10**-6

#Laser Detunings -2*sc.pi*
Deltag = 2*sc.pi*0 #detuning of 493 laser
Deltar = 2*sc.pi*0#detuning of 650 laser
DeltaQ= -2*sc.pi*0.000#detuning of 1762 laser

#Rabi Frequnecies
Om12 = 2*sc.pi*0 #Rabi frequency of |1> to |2> 
Om23 = 2*sc.pi*0 #Rabi frequency of |2> to |3>
Om14 = 2*sc.pi*0.018 #Rabi frequency of |1> to |4>

#Decay rates
gamma21 =  0#2*sc.pi*15.1*MHz #Decay rate of |2> to |1>
gamma23 =  0#2*sc.pi*5.3*MHz #Decay rate of |2> to |3>
gamma41 =  2*sc.pi*10.1*10**-9 #Decay rate of |4> to |1> used a lifetime of 31.2 s

#Laser Linewidths
gammalg = 0#2*sc.pi*2*MHz #493 laser linewidth
gammalr = 0#2*sc.pi*2*MHz #650 laser linewidth
gammalQ = 2*sc.pi*0.002 #1762 laser linewidth


##Initial state of system  
psi0 = basis(4,0)
tlist = np.linspace(0,350,500)

#Operators between |n> and |m> 
sig11 = basis(4,0) * basis(4,0).dag()
sig12 = basis(4,0) * basis(4,1).dag()
sig13 = basis(4,0) * basis(4,2).dag()
sig14 = basis(4,0) * basis(4,3).dag()
sig21 = basis(4,1) * basis(4,0).dag()
sig22 = basis(4,1) * basis(4,1).dag()
sig23 = basis(4,1) * basis(4,2).dag()
sig24 = basis(4,1) * basis(4,3).dag()
sig31 = basis(4,2) * basis(4,0).dag()
sig32 = basis(4,2) * basis(4,1).dag()
sig33 = basis(4,2) * basis(4,2).dag()
sig34 = basis(4,2) * basis(4,3).dag()
sig41 = basis(4,3) * basis(4,0).dag()
sig42 = basis(4,3) * basis(4,1).dag()
sig43 = basis(4,3) * basis(4,2).dag()
sig44 = basis(4,3) * basis(4,3).dag()

##Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
C41 = np.sqrt(gamma41) * sig14 #Decay from |4> to |1>
Clg = np.sqrt(2*gammalg) * sig11 #From 493 laser linewidth
Clr = np.sqrt(2*gammalr) * sig33 #From 650 laser linewidth
ClQ = np.sqrt(2*gammalQ) * sig44 #From 1762 laser linewidth

H0 = (0.5*Om14*sig14 + 0.5*Om14*sig41 + DeltaQ*sig44)

# collapse operators
c_ops_list = [C21,C23,C41,Clg,Clr,ClQ]

output = mesolve(H0, psi0, tlist, c_ops_list, [sig44], {},options=Options(nsteps=5000000))  


fig, ax = plt.subplots(figsize=(8,5))

ax.plot(tlist, output.expect[0], label="D5/2 State Prob")
#ax.plot(tlist, output.expect[1], label="Atom excited state")
ax.legend()
ax.set_xlabel('Time [us]')
ax.set_ylabel('Occupation probability')
ax.set_title('1762 nm Rabi oscillations');

#Save graph data
output_data = np.vstack((tlist, output.expect[0])) # join time and expt˓→data
file_data_store('G:\\Shared drives\\Ions\\Ion Data\\1762 Rabi Data\\June\\fit3.dat', output_data.T, numtype="real") # Note the .T for transpose!

