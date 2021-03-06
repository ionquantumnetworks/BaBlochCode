# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 16:00:16 2017

@author: James
"""
import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from qutip import *

#Units
MHz = 1*10**6
kHz = 1*10**3
ns = 1*10**-9
us = 1*10**-6

#Laser Detunings -2*sc.pi*
Deltag = 2*sc.pi*0*MHz #detuning of 493 laser
Deltar = 2*sc.pi*12*MHz#detuning of 650 laser
DeltaQ= 2*sc.pi*0 #detuning of 1762 laser

#Rabi Frequnecies
Om12 = 2*sc.pi*15*MHz #Rabi frequency of |1> to |2> 
Om23 = 2*sc.pi*5*MHz #Rabi frequency of |2> to |3>
Om14 = 2*sc.pi*0*MHz #Rabi frequency of |1> to |4>

#Decay rates
gamma21 =  2*sc.pi*15.1*MHz #Decay rate of |2> to |1>
gamma23 =  2*sc.pi*5.3*MHz #Decay rate of |2> to |3>
gamma41 =  2*sc.pi*5.1*10**-3 #Decay rate of |4> to |1> used a lifetime of 31.2 s

#Laser Linewidths
gammalg = 0#2*sc.pi*2*MHz #493 laser linewidth
gammalr = 0#2*sc.pi*2*MHz #650 laser linewidth
gammalQ = 2*sc.pi*100*kHz #1762 laser linewidth

tlist = np.linspace(-0*us, 1*us, 1000) #List of points for plotting purposes

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

#operators to input into mesolve
e_opsS = [sig11]
e_opsP = [sig22]
e_opsD = [sig33]
e_opsQ = [sig44]

#Hamiltonian of system RWA
H = (Deltag*sig11 + 0.5*Om12*sig12 + 0.5*Om14*sig14 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32 + Deltar*sig33 + 0.5*Om14*sig41 + DeltaQ*sig44)

##Initial state of system  
psi0 = basis(4,0)

##Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
C41 = np.sqrt(gamma41) * sig14 #Decay from |4> to |1>
Clg = np.sqrt(2*gammalg) * sig11 #From 493 laser linewidth
Clr = np.sqrt(2*gammalr) * sig33 #From 650 laser linewidth
ClQ = np.sqrt(2*gammalQ) * sig44 #From 1762 laser linewidth
c_ops = [C21,C23,C41,Clg,Clr,ClQ]
#
##Solutions
nS = mesolve(H, psi0, tlist, c_ops, e_opsS) #Solves Hamiltionian at point on tlist
nP = mesolve(H, psi0, tlist, c_ops, e_opsP) #Solves Hamiltionian at point on tlist
nD = mesolve(H, psi0, tlist, c_ops, e_opsD) #Solves Hamiltionian at point on tlist
nQ = mesolve(H, psi0, tlist, c_ops, e_opsQ) #Solves Hamiltionian at point on tlist
Spop= np.real(nS.expect[0])
Ppop= np.real(nP.expect[0])
Dpop= np.real(nD.expect[0])
Qpop= np.real(nQ.expect[0])
print('Spop =' + str(Spop[-1]))
print('Ppop =' + str(Ppop[-1]))
print('Dpop =' + str(Dpop[-1]))
print('Qpop =' + str(Qpop[-1]))
#final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
#fexpt = expect(e_ops, final_state) #Calculates expectation values (e_ops list) for Hamiltonion at t = inifinity
#print(H)
fig, ax = plt.subplots(figsize=(12,8))
ax.plot(tlist, Spop, 'b', tlist, Ppop, 'r', tlist, Dpop, 'g', tlist, Qpop, 'y')
ax.set_xlabel('Time')
ax.set_ylabel('Occupation probability')
#ax.set_title('TITLE')
ax.legend(("S state", "P state", "D state", "Qubit D state"), loc=0);
#Save graph data
#output_data = np.vstack((tlist, n.expect)) # join time and expt˓→data
#file_data_store('E:\\IonTrapData\\Tests\\PhotonShape\\8Junephoton.dat', output_data.T, numtype="real") # Note the .T for transpose!

