#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 09:56:49 2018

@author: jmhannegan
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pandas import *


###### Defining Atomic Constants and Laser Parameters ###########
#units MHz -> time units us
#detunings
deltag = -0*(2*np.pi)
deltar = -1 # *(2*np.pi)
detuninglist = [deltag,deltag,0,0,deltar,deltar,deltar,deltar]
#rabi frequencies
Omg = 0#10*(2*np.pi)
Omr = 50*(2*np.pi)
#natural decay rate
gPS = 2*np.pi*15.1
gPD = 2*np.pi*5.3
#laser linewidths
lg = 1*np.pi*2
lr = 1*np.pi*2
#Lande G-Factors
gfactors=[2,2,2/3,2/3,4/5,4/5,4/5,4/5]
mvalues=[-1,1,-1,1,-3,-1,1,3] #mvalues without factor of 2


#Bfield / Zeeman Splitting
uBohr = 1 #Bohr Magneton
B = 5 #Magnetic Field
hbar= 1 #2pi * plank constant
u= uBohr*B

#Branching Ratios



##################################################################

#For now, assume two laser fields, one at 493 nm and the other at 650 nm#
#light is propagating in y-direction#
#The amount of sigma/pi light will be calculated from angle to B-field#
alpha = .2*np.pi #angle alpha

#RWA matrix
#Diagonals
H0 = 0 * basis(8,0) * basis(8,0).dag()
for x in range(0,8):
    H0 = H0 + detuninglist[x]*basis(8,x)*basis(8,x).dag() + 0.5 * mvalues[x]*u*gfactors[x]*basis(8,x)*basis(8,x).dag()
#Off Diagonals
H02= (Omg/np.sqrt(3))*cos(alpha)
H03= -(Omg/np.sqrt(3))*sin(alpha)
H12= H03
H13= H03
H20= H02
H21= H03
H24= -(Omr/2)*sin(alpha)
H25= -(Omr/np.sqrt(3))*cos(alpha)
H26= (Omr/2/np.sqrt(3))*sin(alpha)
H30= H21
H31= H21
H35= -H26
H36= -2*H26
H37= -H24
H42= H24
H52= H25
H53= H35
H62= H26
H63= H36
H73= H37


Hoffdiag= H02*basis(8,0)*basis(8,2).dag() + \
 H03*basis(8,0)*basis(8,3).dag() + \
  H12*basis(8,1)*basis(8,2).dag() + \
   H13*basis(8,1)*basis(8,3).dag() + \
    H20*basis(8,2)*basis(8,0).dag() + \
     H21*basis(8,2)*basis(8,1).dag() + \
      H24*basis(8,2)*basis(8,4).dag() + \
       H25*basis(8,2)*basis(8,5).dag() + \
        H26*basis(8,2)*basis(8,6).dag() + \
         H30*basis(8,3)*basis(8,0).dag() + \
          H31*basis(8,3)*basis(8,1).dag() + \
           H35*basis(8,3)*basis(8,5).dag() + \
            H36*basis(8,3)*basis(8,6).dag() + \
             H37*basis(8,3)*basis(8,7).dag() + \
              H42*basis(8,4)*basis(8,2).dag() + \
               H52*basis(8,5)*basis(8,2).dag() + \
                H53*basis(8,5)*basis(8,3).dag() + \
                 H62*basis(8,6)*basis(8,2).dag() + \
                  H63*basis(8,6)*basis(8,3).dag() + \
                   H73*basis(8,7)*basis(8,3).dag()

Htot=H0+Hoffdiag
#Decay Terms See Oberst Mater's Thesis Innsbruck
C1= np.sqrt(2*gPS/3)*basis(8,0)*basis(8,3).dag()
C2= np.sqrt(2*gPS/3)*basis(8,1)*basis(8,2).dag()
C3= np.sqrt(gPS/3)*(basis(8,0)*basis(8,2).dag()-basis(8,1)*basis(8,3).dag())
C4= np.sqrt(gPD/2)*basis(8,4)*basis(8,2).dag()+np.sqrt(gPD/6)*basis(8,5)*basis(8,3).dag()
C5= np.sqrt(gPD/6)*basis(8,6)*basis(8,2).dag()+np.sqrt(gPD/2)*basis(8,7)*basis(8,3).dag()
C6= np.sqrt(gPD/3)*(basis(8,5)*basis(8,2).dag()-basis(8,6)*basis(8,3).dag())
#linewidth terms
C7= np.sqrt(2*lg)*(basis(8,0)*basis(8,0).dag()+basis(8,1)*basis(8,1).dag())
C8= np.sqrt(2*lr)*(basis(8,4)*basis(8,4).dag()+basis(8,5)*basis(8,5).dag()
    +basis(8,6)*basis(8,6).dag()+basis(8,7)*basis(8,7).dag())

c_ops=[C1,C2,C3,C4,C5,C6,C7,C8]




#####End of Setup#############
tlist = np.linspace(0,.020,1000) #gives times of evaluation (start, stop, # of steps)

s1= basis(8,0)*basis(8,0).dag()
s2= basis(8,1)*basis(8,1).dag()
p1= basis(8,2)*basis(8,2).dag()
p2= basis(8,3)*basis(8,3).dag()
d1= basis(8,4)*basis(8,4).dag()
d2= basis(8,5)*basis(8,5).dag()
d3= basis(8,6)*basis(8,6).dag()
d4= basis(8,7)*basis(8,7).dag()

rho_initial = basis(8,7)*basis(8,7).dag() #Initial state of density matrix
sigSstates= basis(8,0)*basis(8,0).dag() + basis(8,1)*basis(8,1).dag()
sigPstates= basis(8,2)*basis(8,2).dag() + basis(8,3)*basis(8,3).dag()
sigDstates= basis(8,4)*basis(8,4).dag() + basis(8,5)*basis(8,5).dag() + basis(8,6)*basis(8,6).dag() + basis(8,7)*basis(8,7).dag()
pop = mesolve(Htot,rho_initial,tlist,c_ops,[s1,s2,p1,p2,d1,d2,d3,d4]) #solve system for times given by t list

plt.plot(tlist,pop.expect[0], label='S1')
plt.plot(tlist,pop.expect[1], label='S2')
plt.plot(tlist,pop.expect[0]+pop.expect[1], label='Total')
plt.legend()
show()

plt.plot(tlist,pop.expect[2], label='P1')
plt.plot(tlist,pop.expect[3], label='P2')
plt.legend()
show()

plt.plot(tlist,pop.expect[4], label='D1')
plt.plot(tlist,pop.expect[5], label='D2')
plt.plot(tlist,pop.expect[6], label='D3')
plt.plot(tlist,pop.expect[7], label='D4')
plt.legend()
show()






#############New code in progress######################

