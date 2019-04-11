# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 16:00:16 2017

@author: James
"""
import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from qutip import *



#Units of frequnecies are in MZ
#Time units will be in us
Deltag = -2*sc.pi*10 #detuning of 493 laser
Deltar = 2*sc.pi*25 #detuning of 650 laser
Om12 = 2*sc.pi*50 #Rabi frequrency of |1> to |2> 
Om23 = 2*sc.pi*0 #Rabi frequrency of |2> to |3>
gamma21 =  2*sc.pi*15.1 #Decay rate of |2> to |1>
gamma23 =  2*sc.pi*5.3 #Decay rate of |2> to |3>
gammalg = 2*sc.pi*2 #493 laser linewidth
gammalr = 2*sc.pi*2 #650 laser linewidth
tlist = np.linspace(0, 2, 1000) #List of points for plotting purposes

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
e_ops = [sig33] #operators to input into mesolve

#Hamiltonian of system RWA
H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)
#Initial state of system
psi0 = one
#Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
Clg = np.sqrt(2*gammalg) * sig11 #From 493 laser linewidth
Clr = np.sqrt(2*gammalr) * sig33 #From 650 laser linewidth
c_ops = [C21,C23,Clg,Clr]

#Solutions
n = mesolve(H, psi0, tlist, c_ops, e_ops) #Solves Hamiltionian at point on tlist
final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
fexpt = expect(e_ops, final_state) #Calculates expectation values (e_ops list) for Hamiltonion at t = inifinity
#print(sig23)
plot_expectation_values(n, show_legend=True,figsize=(8, 8))

#Save graph data
#output_data = np.vstack((tlist, n.expect)) # join time and expt˓→data
#file_data_store('E:\\IonTrapData\\Tests\\PhotonShape\\8Junephoton.dat', output_data.T, numtype="real") # Note the .T for transpose!

