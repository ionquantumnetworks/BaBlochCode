# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:39:01 2019

@author: John Hannegan
"""

import os
import sys
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from pandas import *
import bariumbloch as bb
import time
qutip.settings.has_mkl=False #I ONLY NEED THIS ON IONS COMPUTER

start=time.time()

#define 8 level basis for now
diags=bb.DiagExpect()

##just pi polarizations along y for now##
lg=0*np.pi*2
lr=0*np.pi*2
alpha=np.pi*50/180#.45
OgPiY=15.1*np.pi*2*(3) #493 Rabi 1
OrPiY=5.3*np.pi*2*(6) #650 Rabi  4
B=2*np.pi*1.4*(20) #B field splitting (in gauss)
Dg= -150*2*np.pi  #50
Dr= -np.pi*2*75 #10


tlist=np.linspace(0,2,1000)
Hoffdiag=bb.BaOffDiag(alpha,0,OgPiY,OrPiY,0,0,0,0,0,0)
c_ops=bb.BaC_ops(lg,lr)

###Vary 650 with 493 Constant#####
begin = -np.pi*250*2
end = np.pi*250*2
diff = end-begin
detuneranger = np.linspace(0, diff, 500)
poplistr=[]
Pstate=bb.DiagExpect()[2]+bb.DiagExpect()[3]

for x in detuneranger:
    H0=bb.BaDiag(Dg,x+begin,B,diags)
    Htot=H0+Hoffdiag
    pop=steadystate(Htot,c_ops)
    poplistr.append(expect(diags[2]+diags[3],pop))
#print(poplist)
plt.plot((detuneranger+begin)/2/np.pi,poplistr,'r')
plt.axis([begin/2/np.pi,end/2/np.pi,0,max(poplistr)*1.1])
show()

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
#    H0=bb.BaDiag(x+begin,Dr,B,diags)
#    Htot=H0+Hoffdiag
#    pop=steadystate(Htot,c_ops)
#    poplist.append(expect(diags[2]+diags[3],pop))
##print(poplist)
#plt.plot((detunerange+begin)/2/np.pi,poplist,'b')
#plt.axis([begin/2/np.pi,end/2/np.pi,0,max(poplist)*1.1])
#show()
#print('650 detuning')
#print(Dr/np.pi/2)
#print('493 rabi')
#print(OgPiY/15.1/np.pi/2)
#print('650 rabi')
#print(OrPiY/5.3/np.pi/2)



#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanpop.dat', 'w')
#for item in poplist:
#  thefile.write("%s\n" % item)
#thefile.close()
#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanfreq.dat', 'w')
#for item in detunerange:
#  thefile.write("%s\n" % (item+begin))
#thefile.close()
#
#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanpopr.dat', 'w')
#for item in poplistr:
#  thefile.write("%s\n" % item)
#thefile.close()
#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\Ba Spectroscopy\\scanfreqr.dat', 'w')
#for item in detuneranger:
#  thefile.write("%s\n" % (item+begin))
#thefile.close()
#
#
#
#
#






#Test to make sure at least one parameter has a nonzero steady state#
#H0=bb.BaDiag(Dg,-np.pi*100,B,diags)
#Htot=H0+Hoffdiag
#pop=mesolve(Htot,diags[1],tlist,c_ops,diags)
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
#print(expect(diags[2]+diags[3],popsteady))
#
#start=time.time()
#H0=bb.BaDiag(Dg,Dr,B,diags)
#print(time.time()-start)
#print('single H0 time')
#
#start=time.time()
#Htot=H0+Hoffdiag
#pop=steadystate(Htot,c_ops)
#print(time.time()-start)
#print('single ss time')