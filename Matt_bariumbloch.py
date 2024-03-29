#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 21:34:31 2019

@author: Matt
"""

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pandas import *
from qutip import basis,steadystate
from pylab import sin,cos

def barium3beam(deltag,deltar,lg,lr,B,alphaY,alphaZ,O493Y,O650Y,O493Z,O650Z,O493sigplus,O493sigminus,O650sigplus,O650sigminus,tlist,rho_initial):
    #natural decay rate
    gPS = 2*np.pi*15.1
    gPD = 2*np.pi*5.3
    #Lande G-Factors
    gfactors=[2,2,2/3,2/3,4/5,4/5,4/5,4/5]
    mvalues=[-1,1,-1,1,-3,-1,1,3] #mvalues without factor of 2
    #Bfield / Zeeman Splitting
    uBohr = 1 #Bohr Magneton
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
    ###########################################################################
    
    ##############Electric Field Definition##############################
    #Magnetic Field in Z direction#
    #One beam for each color along y direction, linear polarization with arbitrary angle to B field given by alphaY#
    #One beam for each color along z direction, with either linear or circular polarization linear field angle given by alphaZ#
    #We will be assuming that each beam of the same color is of the same frequency so that the RWA is easy(/possible?)#
    
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
    ##########
    
    #Total Fields#
    E493tot=Efieldy493+Efieldz493+Efield493sigplus+Efield493sigminus
    E650tot=Efieldy650+Efieldz650+Efield650sigplus+Efield650sigminus
    ##########################################################################
    
    ###################Constructing Hint##################################
    Rabilist=[[16,np.dot(E493tot,D20)],[24,np.dot(E493tot,D30)],[17,np.dot(E493tot,D21)],[25,np.dot(E493tot,D31)],[2,np.conjugate(np.dot(E493tot,D20))],[3,np.conjugate(np.dot(E493tot,D30))],[10,np.conjugate(np.dot(E493tot,D21))],[11,np.conjugate(np.dot(E493tot,D31))],[20,np.dot(E650tot,D24)],[28,np.dot(E650tot,D34)],[21,np.dot(E650tot,D25)],[29,np.dot(E650tot,D35)],[22,np.dot(E650tot,D26)],[30,np.dot(E650tot,D36)],[23,np.dot(E650tot,D27)],[31,np.dot(E650tot,D37)],[34,np.conjugate(np.dot(E650tot,D24))],[35,np.conjugate(np.dot(E650tot,D34))],[42,np.conjugate(np.dot(E650tot,D25))],[43,np.conjugate(np.dot(E650tot,D35))],[50,np.conjugate(np.dot(E650tot,D26))],[51,np.conjugate(np.dot(E650tot,D36))],[58,np.conjugate(np.dot(E650tot,D27))],[59,np.conjugate(np.dot(E650tot,D37))]]
    
    
    Hoffdiaglist=[0]*64
    for x in range(len(Rabilist)):
        Hoffdiaglist[Rabilist[x][0]]=Rabilist[x][1]
    
    #RWA matrix
    #Diagonals
    detuninglist = [deltag,deltag,0,0,deltar,deltar,deltar,deltar]
    H0 = 0 * basis(8,0) * basis(8,0).dag()
    for x in range(0,8):
        H0 = H0 + detuninglist[x]*basis(8,x)*basis(8,x).dag() + 0.5 * mvalues[x]*u*gfactors[x]*basis(8,x)*basis(8,x).dag()
    #offDiags
    Hoffdiag=0*basis(8,0) * basis(8,0).dag()
    for x in range(0,8):
        for y in range(0,8):
            Hoffdiag=Hoffdiag + Hoffdiaglist[8*x+y] * basis(8,x)*basis(8,y).dag()
    #total
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
    C8= np.sqrt(2*lr)*(basis(8,4)*basis(8,4).dag()+basis(8,5)*basis(8,5).dag()\
    +basis(8,6)*basis(8,6).dag()+basis(8,7)*basis(8,7).dag())
    #Putting all into one list
    c_ops=[C1,C2,C3,C4,C5,C6,C7,C8]
    
    s1= basis(8,0)*basis(8,0).dag()
    s2= basis(8,1)*basis(8,1).dag()
    p1= basis(8,2)*basis(8,2).dag()
    p2= basis(8,3)*basis(8,3).dag()
    d1= basis(8,4)*basis(8,4).dag()
    d2= basis(8,5)*basis(8,5).dag()
    d3= basis(8,6)*basis(8,6).dag()
    d4= basis(8,7)*basis(8,7).dag()
    opts = Options(rhs_reuse=True)
    population = mesolve(Htot,rho_initial,tlist,c_ops,[s1,s2,p1,p2,d1,d2,d3,d4],options=opts) #solve system for times given by t list
    #plt.plot(tlist,pop.expect[0], label='S1')
    #plt.plot(tlist,pop.expect[1], label='S2')
    #plt.plot(tlist,pop.expect[0]+pop.expect[1], label='Total')
    #plt.legend()
    #show()
    return population

def barium3beam_ss(deltag,deltar,lg,lr,B,alphaY,alphaZ,O493Y,O650Y,O1762Y,O493Z,O650Z,O1762Z,O493sigplus,O493sigminus,O650sigplus,O650sigminus,O1762sigplus,O1762sigminus):
    #natural decay rate
    gPS = 2*np.pi*15.1
    gPD = 2*np.pi*5.3
    #Lande G-Factors
    gfactors=[2,2,2/3,2/3,4/5,4/5,4/5,4/5]
    mvalues=[-1,1,-1,1,-3,-1,1,3] #mvalues without factor of 2
    #Bfield / Zeeman Splitting
    uBohr = 1 #Bohr Magneton
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
    ###########################################################################
    
    ##############Electric Field Definition##############################
    #Magnetic Field in Z direction#
    #One beam for each color along y direction, linear polarization with arbitrary angle to B field given by alphaY#
    #One beam for each color along z direction, with either linear or circular polarization linear field angle given by alphaZ#
    #We will be assuming that each beam of the same color is of the same frequency so that the RWA is easy(/possible?)#
    
    #y Beams# 
    Efieldy493 = O493Y * np.array([sin(alphaY),0,cos(alphaY)])
    Efieldy650 = O650Y * np.array([sin(alphaY),0,cos(alphaY)])
    Efieldy1762 = O1762Y * np.array([sin(alphaY),0,cos(alphaY)])
    
    #z beams#
    Efieldz493=O493Z * np.array([sin(alphaZ),cos(alphaZ),0])
    Efieldz650=O650Z * np.array([sin(alphaZ),cos(alphaZ),0])
    Efieldz1762=O1762Z * np.array([sin(alphaZ),cos(alphaZ),0])
    
    Efield493sigplus = O493sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield493sigminus = O493sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    Efield650sigplus = O650sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield650sigminus = O650sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    Efield1762sigplus = O1762sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield1762sigminus = O1762sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    ##########
    
    #Total Fields#
    E493tot=Efieldy493+Efieldz493+Efield493sigplus+Efield493sigminus
    E650tot=Efieldy650+Efieldz650+Efield650sigplus+Efield650sigminus
    E1762tot=Efieldy1762+Efieldz1762+Efield1762sigplus+Efield1762sigminus
    ##########################################################################
    
    ###################Constructing Hint##################################
    Rabilist=[[16,np.dot(E493tot,D20)],[24,np.dot(E493tot,D30)],[17,np.dot(E493tot,D21)],[25,np.dot(E493tot,D31)],[2,np.conjugate(np.dot(E493tot,D20))],[3,np.conjugate(np.dot(E493tot,D30))],[10,np.conjugate(np.dot(E493tot,D21))],[11,np.conjugate(np.dot(E493tot,D31))],[20,np.dot(E650tot,D24)],[28,np.dot(E650tot,D34)],[21,np.dot(E650tot,D25)],[29,np.dot(E650tot,D35)],[22,np.dot(E650tot,D26)],[30,np.dot(E650tot,D36)],[23,np.dot(E650tot,D27)],[31,np.dot(E650tot,D37)],[34,np.conjugate(np.dot(E650tot,D24))],[35,np.conjugate(np.dot(E650tot,D34))],[42,np.conjugate(np.dot(E650tot,D25))],[43,np.conjugate(np.dot(E650tot,D35))],[50,np.conjugate(np.dot(E650tot,D26))],[51,np.conjugate(np.dot(E650tot,D36))],[58,np.conjugate(np.dot(E650tot,D27))],[59,np.conjugate(np.dot(E650tot,D37))]]
    
    
    Hoffdiaglist=[0]*64
    for x in range(len(Rabilist)):
        Hoffdiaglist[Rabilist[x][0]]=Rabilist[x][1]
    
    #RWA matrix
    #Diagonals
    detuninglist = [deltag,deltag,0,0,deltar,deltar,deltar,deltar]
    H0 = 0 * basis(8,0) * basis(8,0).dag()
    for x in range(0,8):
        H0 = H0 + detuninglist[x]*basis(8,x)*basis(8,x).dag() + 0.5 * mvalues[x]*u*gfactors[x]*basis(8,x)*basis(8,x).dag()
    #offDiags
    Hoffdiag=0*basis(8,0) * basis(8,0).dag()
    for x in range(0,8):
        for y in range(0,8):
            Hoffdiag=Hoffdiag + Hoffdiaglist[8*x+y] * basis(8,x)*basis(8,y).dag()
    #total
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
    C8= np.sqrt(2*lr)*(basis(8,4)*basis(8,4).dag()+basis(8,5)*basis(8,5).dag()\
    +basis(8,6)*basis(8,6).dag()+basis(8,7)*basis(8,7).dag())
    #Putting all into one list
    c_ops=[C1,C2,C3,C4,C5,C6,C7,C8]
    opts = Options(rhs_reuse=True)
    
    population_ss = steadystate(Htot,c_ops) #solve system for times given by t list
    #plt.plot(tlist,pop.expect[0], label='S1')
    #plt.plot(tlist,pop.expect[1], label='S2')
    #plt.plot(tlist,pop.expect[0]+pop.expect[1], label='Total')
    #plt.legend()
    #show()
    return population_ss

    
    
def BaOffDiag(alphaY,alphaZ,O493Y,O650Y,O1762Y,O493Z,O650Z,O1762Z,O493sigplus,O493sigminus,O650sigplus,O650sigminus,O1762sigplus, O1762sigminus):
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
    ###########################################################################
    
    ##############Electric Field Definition##############################
    #Magnetic Field in Z direction#
    #One beam for each color along y direction, linear polarization with arbitrary angle to B field given by alphaY#
    #One beam for each color along z direction, with either linear or circular polarization linear field angle given by alphaZ#
    #We will be assuming that each beam of the same color is of the same frequency so that the RWA is easy(/possible?)#
    
    #y Beams#
    Efieldy493 = O493Y * np.array([sin(alphaY),0,cos(alphaY)])
    Efieldy650 = O650Y * np.array([sin(alphaY),0,cos(alphaY)])
    Efieldy1762 = O1762Y * np.array([sin(alphaY),0,cos(alphaY)])
                                  
    #z beams#
    Efieldz493=O493Z * np.array([sin(alphaZ),cos(alphaZ),0])
    Efieldz650=O650Z * np.array([sin(alphaZ),cos(alphaZ),0])
    Efieldz1762=O1762Z * np.array([sin(alphaZ),cos(alphaZ),0])
    
    #Circular Beams#
    Efield493sigplus = O493sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield493sigminus = O493sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    Efield650sigplus = O650sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield650sigminus = O650sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    Efield1762sigplus = O1762sigplus/np.sqrt(2) * np.array([1, 1j, 0])
    Efield1762sigminus = O1762sigminus/np.sqrt(2) * np.array([1, -1j, 0])
    ##########
    
    #Total Fields#
    E493tot=Efieldy493+Efieldz493+Efield493sigplus+Efield493sigminus
    E650tot=Efieldy650+Efieldz650+Efield650sigplus+Efield650sigminus
    E1762tot=Efieldy1762+Efieldz1762+Efield1762sigplus+Efield1762sigminus
    ##########################################################################
    
    ###################Constructing Hint##################################
    Rabilist=[[16,np.dot(E493tot,D20)],[24,np.dot(E493tot,D30)],[17,np.dot(E493tot,D21)],[25,np.dot(E493tot,D31)],[2,np.conjugate(np.dot(E493tot,D20))],[3,np.conjugate(np.dot(E493tot,D30))],[10,np.conjugate(np.dot(E493tot,D21))],[11,np.conjugate(np.dot(E493tot,D31))],[20,np.dot(E650tot,D24)],[28,np.dot(E650tot,D34)],[21,np.dot(E650tot,D25)],[29,np.dot(E650tot,D35)],[22,np.dot(E650tot,D26)],[30,np.dot(E650tot,D36)],[23,np.dot(E650tot,D27)],[31,np.dot(E650tot,D37)],[34,np.conjugate(np.dot(E650tot,D24))],[35,np.conjugate(np.dot(E650tot,D34))],[42,np.conjugate(np.dot(E650tot,D25))],[43,np.conjugate(np.dot(E650tot,D35))],[50,np.conjugate(np.dot(E650tot,D26))],[51,np.conjugate(np.dot(E650tot,D36))],[58,np.conjugate(np.dot(E650tot,D27))],[59,np.conjugate(np.dot(E650tot,D37))]]
    
    
    Hoffdiaglist=[0]*64
    for x in range(len(Rabilist)):
        Hoffdiaglist[Rabilist[x][0]]=Rabilist[x][1]

    #offDiagsRWA
    Hoffdiag=0*basis(8,0) * basis(8,0).dag()
    for x in range(0,8):
        for y in range(0,8):
            Hoffdiag=Hoffdiag + Hoffdiaglist[8*x+y] * basis(8,x)*basis(8,y).dag()

    return Hoffdiag
    
def DiagExpect():
    s1= basis(8,0)*basis(8,0).dag()
    s2= basis(8,1)*basis(8,1).dag()
    p1= basis(8,2)*basis(8,2).dag()
    p2= basis(8,3)*basis(8,3).dag()
    d1= basis(8,4)*basis(8,4).dag()
    d2= basis(8,5)*basis(8,5).dag()
    d3= basis(8,6)*basis(8,6).dag()
    d4= basis(8,7)*basis(8,7).dag()
    return [s1,s2,p1,p2,d1,d2,d3,d4]

def BaDiag(deltag,deltar,B,basislist):
    #Lande G-Factors
    gfactors=[2,2,2/3,2/3,4/5,4/5,4/5,4/5]
    mvalues=[-1,1,-1,1,-3,-1,1,3] #mvalues without factor of 2
    #Bfield / Zeeman Splitting
    uBohr = 1 #Bohr Magneton
    u= uBohr*B
    #RWA matrix
    #Diagonals
    detuninglist = [deltag,deltag,0,0,deltar,deltar,deltar,deltar]
    H0 = 0
    for x in range(0,8):
        H0 = H0 + (detuninglist[x] + 0.5 * mvalues[x]*u*gfactors[x])*basislist[x]
    return H0
    
def BaC_ops(lg,lr,ll):
    #natural decay rate
    gPS = 2*np.pi*15.1
    gPD = 2*np.pi*5.3
    #Decay Terms See Oberst Mater's Thesis Innsbruck
    C1= np.sqrt(2*gPS/3)*basis(8,0)*basis(8,3).dag()
    C2= np.sqrt(2*gPS/3)*basis(8,1)*basis(8,2).dag()
    C3= np.sqrt(gPS/3)*(basis(8,0)*basis(8,2).dag()-basis(8,1)*basis(8,3).dag())
    C4= np.sqrt(gPD/2)*basis(8,4)*basis(8,2).dag()+np.sqrt(gPD/6)*basis(8,5)*basis(8,3).dag()
    C5= np.sqrt(gPD/6)*basis(8,6)*basis(8,2).dag()+np.sqrt(gPD/2)*basis(8,7)*basis(8,3).dag()
    C6= np.sqrt(gPD/3)*(basis(8,5)*basis(8,2).dag()-basis(8,6)*basis(8,3).dag())
    #linewidth terms
    C7= np.sqrt(2*lg)*(basis(8,0)*basis(8,0).dag()+basis(8,1)*basis(8,1).dag())
    C8= np.sqrt(2*lr)*(basis(8,4)*basis(8,4).dag()+basis(8,5)*basis(8,5).dag()\
    +basis(8,6)*basis(8,6).dag()+basis(8,7)*basis(8,7).dag())
    #Putting all into one list
    c_ops=[C1,C2,C3,C4,C5,C6,C7,C8]
    
    return c_ops
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#tlist=np.linspace(0,20,1000)
#initialstate=basis(8,5)*basis(8,5).dag()

#pop=barium3beam(2*np.pi*10,-2*np.pi*10,0,0,10,np.pi/3,0,30*np.pi,30*np.pi,0,0,0,0,0,0,timelist,initialstate)
#print(pop)

#plt.plot(tlist,pop.expect[0], label='S1')
#plt.plot(tlist,pop.expect[1], label='S2')
#plt.plot(tlist,pop.expect[0]+pop.expect[1], label='Total')
#plt.legend()
#show()

#plt.plot(tlist,pop.expect[2], label='P1')
#plt.plot(tlist,pop.expect[3], label='P2')
#plt.plot(tlist,pop.expect[2]+pop.expect[3],label='tot')
#plt.legend()
#show()

#plt.plot(tlist,pop.expect[4], label='D1')
#plt.plot(tlist,pop.expect[5], label='D2')
#plt.plot(tlist,pop.expect[6], label='D3')
#plt.plot(tlist,pop.expect[7], label='D4')
#plt.plot(tlist,pop.expect[4]+pop.expect[5]+pop.expect[6]+pop.expect[7],label='tot')
#plt.legend()
#show()

#pop_ss=barium3beam_ss(2*np.pi*10,-2*np.pi*10,0,0,10,np.pi/3,0,30*np.pi,30*np.pi,0,0,0,0,0,0)
#print(expect(basis(8,1)*basis(8,1).dag(),pop_ss))
#print(pop.expect[1][-1])