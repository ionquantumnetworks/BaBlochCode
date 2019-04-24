# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:53:52 2019

@author: John Hannegan
"""


from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pandas import *

###### Defining Atomic Constants and other Laser Parameters ###########
#units MHz -> time units us
#detunings
deltag = -20*(2*np.pi)
deltar = 0 *(2*np.pi)
detuninglist = [deltag,deltag,0,0,deltar,deltar,deltar,deltar]
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
B = 10 #Magnetic Field
hbar= 1 #2pi * plank constant
u= uBohr*B


###############Barium Dipole Matrix Operators############################

##S1/2 to P1/2##
D20=np.array([0,0,1/np.sqrt(3)])
D30=np.array([-1/np.sqrt(3),1j/np.sqrt(3),0])
D21=np.array([-1/np.sqrt(3),-1j/np.sqrt(3),0])
D31=np.array([0,0,-1/np.sqrt(3)])
###############

##D3/2 to P1/2##
D24=np.array([-1/2,1j/2,0])
D34=np.array([0,0,0])
D25=np.array([0,0,-1/np.sqrt(3)])
D35=np.array([-1/2/np.sqrt(3),1j/2/np.sqrt(3),0])
D26=np.array([1/2/np.sqrt(3),1j/2/np.sqrt(3),0])
D36=np.array([0,0,-1/np.sqrt(3)])
D27=np.array([0,0,0])
D37=np.array([1/2,1j/2,0])

########################

##############Electric Field Definition##############################
#Magnetic Field in Z direction#
#One beam for each color along y direction, linear polarization with arbitrary angle to B field given by alphaY#
#One beam for each color along z direction, with either linear or circular polarization linear field angle given by alphaZ#
#We will be assuming that each beam of the same color is of the same frequency so that the RWA is easy(/possible?)#

#Polarization angles of Pi beams#
alphaY = np.pi/2
alphaZ = np.pi*0

#Pi Rabifreq/efield for each Direction#
O493Y=2*np.pi*0
O493Z=2*np.pi*10

O650Y=2*np.pi*0
O650Z=2*np.pi*10

#Sigma Rabifreq/efields along z#
O493sigplus = 2*np.pi*0
O493sigminus = 2*np.pi*0

O650sigplus = 0
O650sigminus = 0


#y Beams#
Efieldy493 = O493Y * np.array([sin(alphaY),0,cos(alphaY)])
Efieldy650 = O650Y * np.array([sin(alphaY),0,cos(alphaY)])

#z beams#
Efieldz493=O493Z * np.array([sin(alphaZ),cos(alphaZ),0])
Efieldz650=O650Z * np.array([sin(alphaZ),cos(alphaZ),0])

Efield493sigplus = O493sigplus/np.sqrt(2) * np.array([1, 1j, 0])
Efield493sigminus = O493sigminus/np.sqrt(2) * np.array([1, -1j, 0])
Efield650sigplus = O650sigplus/np.sqrt(2) * np.array([1, 1j, 0])
Efield650sigminus = O650sigminus/np.sqrt(2) * np.array([1, -1j, 0])

#Total Fields#
E493tot=Efieldy493+Efieldz493+Efield493sigplus+Efield493sigminus
E650tot=Efieldy650+Efieldz650+Efield650sigplus+Efield650sigminus
######################################################################

###################Constructing Hint##################################
Rabilist=[[16,np.dot(E493tot,D20)],[24,np.dot(E493tot,D30)],[17,np.dot(E493tot,D21)],[25,np.dot(E493tot,D31)],[2,np.conjugate(np.dot(E493tot,D20))],[3,np.conjugate(np.dot(E493tot,D30))],[10,np.conjugate(np.dot(E493tot,D21))],[11,np.conjugate(np.dot(E493tot,D31))],[20,np.dot(E650tot,D24)],[28,np.dot(E650tot,D34)],[21,np.dot(E650tot,D25)],[29,np.dot(E650tot,D35)],[22,np.dot(E650tot,D26)],[30,np.dot(E650tot,D36)],[23,np.dot(E650tot,D27)],[31,np.dot(E650tot,D37)],[34,np.conjugate(np.dot(E650tot,D24))],[35,np.conjugate(np.dot(E650tot,D34))],[42,np.conjugate(np.dot(E650tot,D25))],[43,np.conjugate(np.dot(E650tot,D35))],[50,np.conjugate(np.dot(E650tot,D26))],[51,np.conjugate(np.dot(E650tot,D36))],[58,np.conjugate(np.dot(E650tot,D27))],[59,np.conjugate(np.dot(E650tot,D37))]]

Hoffdiaglist=[0]*64

for x in range(len(Rabilist)):
    Hoffdiaglist[Rabilist[x][0]]=Rabilist[x][1]

#print(Hoffdiaglist)


#Htest=[0]*64
#Htest[1]=Rabilist[1][1]
#print(Htest)
#print(np.dot(Efieldz493,D30))
#print(np.dot([1,1],[1,1j]))
#print(E493tot)
#print(sin(alphaY))

##################################################################

#RWA matrix
#Diagonals
H0 = 0 * basis(8,0) * basis(8,0).dag()
for x in range(0,8):
    H0 = H0 + detuninglist[x]*basis(8,x)*basis(8,x).dag() + 0.5 * mvalues[x]*u*gfactors[x]*basis(8,x)*basis(8,x).dag()

Hoffdiag=0*basis(8,0) * basis(8,0).dag()
for x in range(0,8):
    for y in range(0,8):
        Hoffdiag=Hoffdiag + Hoffdiaglist[8*x+y] * basis(8,x)*basis(8,y).dag()
        #print(x,y)

#print(Hoffdiaglist[8*7+3])
#print(Hoffdiag[7,3])


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
tlist = np.linspace(0,20,1000) #gives times of evaluation (start, stop, # of steps)

s1= basis(8,0)*basis(8,0).dag()
s2= basis(8,1)*basis(8,1).dag()
p1= basis(8,2)*basis(8,2).dag()
p2= basis(8,3)*basis(8,3).dag()
d1= basis(8,4)*basis(8,4).dag()
d2= basis(8,5)*basis(8,5).dag()
d3= basis(8,6)*basis(8,6).dag()
d4= basis(8,7)*basis(8,7).dag()

rho_initial = basis(8,5)*basis(8,5).dag() #Initial state of density matrix
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
plt.plot(tlist,pop.expect[2]+pop.expect[3],label='tot')
plt.legend()
show()

plt.plot(tlist,pop.expect[4], label='D1')
plt.plot(tlist,pop.expect[5], label='D2')
plt.plot(tlist,pop.expect[6], label='D3')
plt.plot(tlist,pop.expect[7], label='D4')
plt.plot(tlist,pop.expect[4]+pop.expect[5]+pop.expect[6]+pop.expect[7],label='tot')
plt.legend()
show()

print(pop)

#plt.plot(tlist,pop.expect[1]+pop.expect[2]+pop.expect[3]+pop.expect[4]+pop.expect[5]+pop.expect[6]+pop.expect[7]+pop.expect[0],label='total pop')
#plt.legend()
#plt.axis([0,.02,0,2])
#show()