# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:08:08 2019

@author: John Hannegan
"""

#this code is meant to be simple three level code to distinguish individual photon emissions when 
#producing a single photon from a 3 level ion with a strong branching ratio

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *


#Turn on or off coupling between separate ions 0=off 1=on
coupling=1

###### Defining Atomic Constants and Laser Parameters ###########
#units MHz -> time units us
#detunings
d12 = 0*(2*np.pi)
d23 = 0*(2*np.pi)
#rabi frequencies
Om12 = 0*(2*np.pi)
Om23 = 5.3*(2*np.pi)*(0)
#natural decay rate
g12 = 2*np.pi*15.1
g23 = 2*np.pi*5.3
#laser linewidths
l12 = 0*np.pi*2
l23 = 0*np.pi*2
##################################################################

######################Setting up Hamiltonian######################
#relevant operators - Two Three Level Systems#

#first ion internal couplings
sig11 = basis(6,0) * basis(6,0).dag()
sig22 = basis(6,1) * basis(6,1).dag()
sig33 = basis(6,2) * basis(6,2).dag()
sig12 = basis(6,0) * basis(6,1).dag()
sig21 = basis(6,1) * basis(6,0).dag()
sig13 = basis(6,0) * basis(6,2).dag()
sig31 = basis(6,2) * basis(6,0).dag()
sig23 = basis(6,1) * basis(6,2).dag()
sig32 = basis(6,2) * basis(6,1).dag()

#second ion internal couplings
sig44 = basis(6,3) * basis(6,3).dag()
sig55 = basis(6,4) * basis(6,4).dag()
sig66 = basis(6,5) * basis(6,5).dag()
sig45 = basis(6,3) * basis(6,4).dag()
sig54 = basis(6,4) * basis(6,3).dag()
sig46 = basis(6,3) * basis(6,5).dag()
sig64 = basis(6,5) * basis(6,3).dag()
sig56 = basis(6,4) * basis(6,5).dag()
sig65 = basis(6,5) * basis(6,4).dag()

#first ion to second ion coupling#
sig62=basis(6,5) * basis(6,1).dag()

###########Incoherent Loss Terms/Couping to environement###########
#spontaneous decay
#within ion 1
C12=np.sqrt(g12)*sig12
C23=(coupling-1)*np.sqrt(g23)*sig32

#within ion 2
C45=np.sqrt(g12)*sig45
C56=(1)*np.sqrt(g23)*sig65

#between ions
C26=(coupling)*np.sqrt(g23)*sig62

#laser linewidth contribution
#ion 1
C11=np.sqrt(l12)*sig11
C33=np.sqrt(l23)*sig33
#ion 2
C44=np.sqrt(l12)*sig44
C66=np.sqrt(l23)*sig66


#putting together collapse operators
c_ops=[C12,C23,C11,C33,C45,C56,C26,C44,C66]
#################################################################


#Hamiltonian
H = d12*sig11 + Om12/2*(sig12+sig21)+ d23*sig33 + Om23/2*(sig23+sig32) \
    + d12*sig44 + Om12/2*(sig45+sig54)+ d23*sig66 + Om23/2*(sig56+sig65)

######## evaluate populations as a function of time ############
tmax = 1 #microseconds
tstep = 10000
tlist = np.linspace(0,tmax,tstep) #gives times of evaluation (start, stop, # of steps)

rho_initial = sig11 * 0 + sig22 * 1 + sig33 * 0 #Initial state of density matrix

pop = mesolve(H,rho_initial,tlist,c_ops,[sig11,sig22,sig33,sig44,sig55,sig66]) #solve system for times given by t list

#plt.plot(tlist,pop.expect[0], label='S State 1')
plt.plot(tlist,pop.expect[1], label='P State 1')
#plt.plot(tlist,pop.expect[2], label='D State 1')
#plt.plot(tlist,pop.expect[3], label='S State 2')
plt.plot(tlist,pop.expect[4]*1, label='P State 2')
plt.plot(tlist,pop.expect[4]+pop.expect[1], label='P State Sum')
#plt.plot(tlist,pop.expect[5], label='D State 2')
plt.legend()
plt.xlim(0,.04)
show()

#plt.plot(tlist,pop.expect[0], label='S State 1')
#plt.plot(tlist,pop.expect[1], label='P State 1')
#plt.plot(tlist,pop.expect[2], label='D State 1')
#plt.plot(tlist,pop.expect[3], label='S State 2')
#plt.plot(tlist,pop.expect[4], label='P State 2')
#plt.plot(tlist,pop.expect[0]+pop.expect[3], label='S State Sum')
#plt.plot(tlist,pop.expect[5], label='D State 2')
#plt.legend()
#show()


plt.plot(np.linspace(0,tstep/tmax/2/np.pi,tstep),np.absolute(np.fft.fft(pop.expect[1]))/np.absolute(np.fft.fft(pop.expect[1]))[0])
plt.xlim(0,50)
show()
plt.plot(np.fft.fftfreq(len(pop.expect[1]))*tstep/2/np.pi,np.absolute(np.fft.fft(pop.expect[1]))/np.absolute(np.fft.fft(pop.expect[1]))[0])
plt.xlim(-0,50)
show()
#print(len(pop.expect[1]))
#print(np.fft.fft((pop.expect[1])))
#print(np.linspace(0,tstep/tmax))