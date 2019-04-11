# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 13:35:12 2017

@author: James
"""
import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from qutip import *
from pylab import *


#Units of frequnecies are in MZ
#Time units will be in us
Deltag = -2*sc.pi*30 #detuning of 493 laser
Deltar = 2*sc.pi*10 #detuning of 650 laser
Om12 = 2*sc.pi*47 #Rabi frequrency of |1> to |2> 
Om23 = 2*sc.pi*30 #Rabi frequrency of |2> to |3>
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

#Hamiltonian of system RWA
H = (Deltag*sig11 + Deltar*sig33 + 0.5*Om12*sig12 + 0.5*Om12*sig21 + 0.5*Om23*sig23 + 0.5*Om23*sig32)

#Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
c_ops = [C21,C23]
e_ops = [sig11,sig22,sig33] #operators to input into mesolve

# initialise plot and line
#x = []
#y = []
tlist = np.linspace(0, 0.2,400) #List of points for plotting purposes

rho0 = one
# first calculate the occupation number as a function of time
prob = mesolve(H, rho0, tlist, c_ops, sig22).expect[0]
# calculate the occupation of the excited state at t = oo
final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
expt = expect(sig22, final_state) #Calculates expectation values for Hamiltonion at t = inifinity in excited state
g = [a / expt for a in prob]
#plot_expectation_values(n, show_legend=True,figsize=(8, 8))
#print(g)

# Plot the solution
plot((tlist[0:40]*1000), g[0:40], 'b',
     (tlist[0:40]*-1000), g[0:40], 'b')
xlabel('Delay Time, $\tau$ [ns]')
ylabel('$g^{(2)}(\tau)$')
show()

plot((tlist[0:400]*1000), g[0:400], 'b',
     (tlist[0:400]*-1000), g[0:400], 'b')
xlabel('Delay Time, $\tau$ [ns]')
ylabel('$g^{(2)}(\tau)$')
show()

#john edit for g1 and spectrum
tlist2=np.linspace(0,1,10000)
corr=correlation_ss(H,tlist2,c_ops,sig21,sig12)
wlist,spec=spectrum_correlation_fft(tlist2,corr)

#other way
#wlist2=np.linspace(0,1000,100000)*2*sc.pi
#spec2=spectrum(H,wlist2,c_ops,sig12.dag(),sig12)
#print(corr)
#plot
#plot(tlist2,(corr/corr[0]).real)
#plot(tlist2,(corr/corr[0]).imag)
#show()
#print()
plot(wlist[0:100],spec[0:100]/spec[0])
show()
plot(wlist,spec/spec[0])
show()
#plot(wlist2,spec2)
#plot(wlist,spec)
#fig, ax = plot
#ax.plot(wlist/(2*sc.pi), spec1, 'b', lw=2, label='me+fft method')
#ax.legend()
#ax.set_xlabel('frequency')
#ax.set_ylabel(' Power Spectrum')
#ax.set_title('Ion Sperctrum')
#ax.set_xlim(wlist[0]/(2*sc.pi),wlist[-1]/(2*sc.pi))
#plt.show()

#thefile = open('G:\\Team Drives\\Ions\\Ba-Yb Blade Trap\\3_Lab Books\\g2\\QutipData\\Ba3lvlg2.dat', 'w')
#for item in g:
#  thefile.write("%s\n" % item)
#thefile.close()