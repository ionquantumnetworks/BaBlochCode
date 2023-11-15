# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:39:01 2019

@author: Matt
"""

import os
import sys
from qutip import *
from qutip import expect, steadystate
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pandas import *
import Matt_bariumbloch as bb
import time
qutip.settings.has_mkl=False #I ONLY NEED THIS ON IONS COMPUTER

start=time.time()

#define 8 level basis for now
HDiags=bb.DiagExpect() #Builds the diagonal hamiltonain basis. A unit matrix

DeltaBlue= 2*np.pi*(-10)  #493 Laser detuning
DeltaRed= -np.pi*2*75  #650 Laser detuning
DeltaLong= (2*np.pi)*5 #1762 Laser detuning

BluePiGamma_Z=(15.1/4)*(np.pi*2) #Rabi frequrency of |1> to |3> and |2> to |4> from Z
BluePiGamma_Y=(15.1/4)*(np.pi*2) #Rabi frequrency of |1> to |3> and |2> to |4> from Y
BlueSigmaPlusGamma=(15.1/4)*(2*np.pi) #Rabi frequrency of |1> to |4>
BlueSigmaMinusGamma=(15.1/4)*(2*np.pi) #Rabi frequrency of |2> to |3>

RedPiGamma_Z=5.3*(np.pi*2) #Rabi frequrency of |6> to |3> and |7> to |4> from Z
RedPiGamma_Y=0*5.3*(np.pi*2) #Rabi frequrency of |6> to |3> and |7> to |4> from Y
RedSigmaPlusGamma=5.3*(2*np.pi) #Rabi frequrency of |5> to |3> and |6> to |4>
RedSigmaMinusGamma=5.3*(2*np.pi) #Rabi frequrency of |7> to |3> and |8> to |4>

LongPiGamma_Z=(15.1/4)*(np.pi*2) #Rabi frequrency of |2> to |7> and |1> to |6> from Z
LongPiGamma_Y=0*(15.1/4)*(np.pi*2) #Rabi frequrency of |2> to |7> and |1> to |6> from Y
LongSigmaPlusGamma=(15.1/4)*(2*np.pi) #Rabi frequrency of |1> to |7>
LongSigmaMinusGamma=(15.1/4)*(2*np.pi) #Rabi frequrency of |2> to |6>

B=2*np.pi*5 #B field splitting (in gauss)

GammaLaserBlue=0*(np.pi*2) #493 Laser Line Width
GammaLaserRed=0*(np.pi*2) #650 Laser Line Width
GammaLaserLong=0*(2*np.pi) #1762 Laser Line Width

alphaY=np.pi*30/180 #Angle from Y
alphaZ=np.pi*30/180 #Angle from Z




tlist=np.linspace(0,2,1000)  #0 to 2 with 1000 elements
HOffDiag=bb.BaOffDiag(alphaY,alphaZ,BluePiGamma_Y,RedPiGamma_Y,LongPiGamma_Y,BluePiGamma_Z,RedPiGamma_Z,LongPiGamma_Z,BlueSigmaPlusGamma,BlueSigmaMinusGamma,RedSigmaPlusGamma,RedSigmaMinusGamma,LongSigmaPlusGamma, LongSigmaMinusGamma) #Generates off diagonal hamiltonian elements. 
c_ops=bb.BaC_ops(GammaLaserBlue,GammaLaserRed,GammaLaserLong) #Makes the C decay operators based on laser linewidth. (The light gamma's are defined in BaC_ops not here)

###Vary 650 with 493 Constant#####
begin = -2*np.pi*60 #Starting Range Of Plot
end = 2*np.pi*50 #Ending Range Of Plot
diff = end-begin 
detuneranger = np.linspace(0, diff, 500) #0 to difference with 500 elements
poplistr=[] #empty variable
Pstate=bb.DiagExpect()[2]+bb.DiagExpect()[3] #The P states basis. Remember Python counts from 0 so [2] and [3] mean
#the |3> and |4> states 

for x in detuneranger: #A loop to go through detuneranger
    H0=bb.BaDiag(DeltaBlue,x+begin,B,HDiags)  #Creates Diagonal Hamiltonian Elements
    Htot=H0+HOffDiag #Total Hamiltonian
    pop=steadystate(Htot,c_ops) #Finds density matrix at t=infinity
    poplistr.append(expect(HDiags[2]+HDiags[3],pop)) #Adds value to poplistr, then repeats loop
#print(poplist)
plt.plot((detuneranger+begin)/2/np.pi,poplistr,'r')
plt.axis([begin/2/np.pi,end/2/np.pi,0,max(poplistr)*1.1])
plt.xlabel('ΔRed [MHz]'); plt.ylabel('Population in |3>+|4>'); plt.title('Hi'); plt.show()



print(time.time()-start)

####Vary 493 with 650 Constant#####
#begin = -np.pi*250*2
#end = np.pi*10*2
#diff = end-begin
#detunerange = np.linspace(0, diff, 500)
#poplist=[]
##Pstate=bb.DiagExpect()[2]+bb.DiagExpect()[3]
#
#for x in detunerange:
#    H0=bb.BaDiag(x+begin,DeltaRed,B,HDiags)
#    Htot=H0+HOffDiag
#    pop=steadystate(Htot,c_ops)
#    poplist.append(expect(HDiags[2]+HDiags[3],pop))
##print(poplist)
#plt.plot((detunerange+begin)/2/np.pi,poplist,'b')
#plt.axis([begin/2/np.pi,end/2/np.pi,0,max(poplist)*1.1])
#show()
#print('650 detuning')
#print(DeltaRed/np.pi/2)
#print('493 rabi')
#print(BluePiGamma_Y/15.1/np.pi/2)
#print('650 rabi')
#print(RedPiGamma_Y/5.3/np.pi/2)



#thefile = open('G:\\Team DeltaRedives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanpop.dat', 'w')
#for item in poplist:
#  thefile.write("%s\n" % item)
#thefile.close()
#thefile = open('G:\\Team DeltaRedives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanfreq.dat', 'w')
#for item in detunerange:
#  thefile.write("%s\n" % (item+begin))
#thefile.close()
#
#thefile = open('G:\\Team DeltaRedives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanpopr.dat', 'w')
#for item in poplistr:
#  thefile.write("%s\n" % item)
#thefile.close()
#thefile = open('G:\\Team DeltaRedives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanfreqr.dat', 'w')
#for item in detuneranger:
#  thefile.write("%s\n" % (item+begin))
#thefile.close()
#
#
#
#
#






#Test to make sure at least one parameter has a nonzero steady state#
#H0=bb.BaDiag(DeltaBlue,-np.pi*100,B,HDiags)
#Htot=H0+HOffDiag
#pop=mesolve(Htot,HDiags[1],tlist,c_ops,HDiags)
#popsteady=steadystate(Htot,c_ops)
#plt.plot(tlist,pop.expect[0], label='S1')
#plt.plot(tlist,pop.expect[1], label='S2')
#plt.plot(tlist,pop.expect[0]+pop.expect[1], label='Total')
#plt.legend()
#show()
#
#plt.plot(tlist,pop.expect[2], label='P1')
#plt.plot(tlist,pop.expect[3], label='P2')
#plt.plot(tlist,pop.expect[2]+pop.expect[3],label='tot')
#plt.legend()
#show()
#
#plt.plot(tlist,pop.expect[4], label='D1')
#plt.plot(tlist,pop.expect[5], label='D2')
#plt.plot(tlist,pop.expect[6], label='D3')
#plt.plot(tlist,pop.expect[7], label='D4')
#plt.plot(tlist,pop.expect[4]+pop.expect[5]+pop.expect[6]+pop.expect[7],label='tot')
#plt.legend()
#show()
#print(expect(HDiags[2]+HDiags[3],popsteady))
#
#start=time.time()
#H0=bb.BaDiag(DeltaBlue,DeltaRed,B,HDiags)
#print(time.time()-start)
#print('single H0 time')
#
#start=time.time()
#Htot=H0+HOffDiag
#pop=steadystate(Htot,c_ops)
#print(time.time()-start)
#print('single ss time')