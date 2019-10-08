# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 18:48:05 2019

@author: James
"""
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
import time
import scipy.constants as sc

def qubit_integrate(Om14, DeltaQ, A,C41,ClQ, psi0, tlist):

    # Hamiltonian
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

    H0 = (0.5*Om14*sig14 + 0.5*Om14*sig41 + DeltaQ*sig44)
    H1 =  A*sig44 #A is sweep rate

    # collapse operators
    c_ops_list = [C21,C23,C41,Clg,Clr,ClQ]

    # evolve and calculate expectation values
    # method 2: a function callback that returns the coefficient for a qobj
    H = [H0, [H1, lambda x,y: x]]
    output = mesolve(H, psi0, tlist, c_ops_list, [sig44], {})  
 

    return output.expect[0]

#Units
MHz = 1*10**6
kHz = 1*10**3
ns = 1*10**-9
us = 1*10**-6

#Laser Detunings -2*sc.pi*
Deltag = 2*sc.pi*0 #detuning of 493 laser
Deltar = 2*sc.pi*0#detuning of 650 laser
DeltaQ= -2*sc.pi*1 #detuning of 1762 laser

#Rabi Frequnecies
Om12 = 2*sc.pi*0 #Rabi frequency of |1> to |2> 
Om23 = 2*sc.pi*0 #Rabi frequency of |2> to |3>
Om14 = 2*sc.pi*1.043 #Rabi frequency of |1> to |4>

#Decay rates
gamma21 =  0#2*sc.pi*15.1*MHz #Decay rate of |2> to |1>
gamma23 =  0#2*sc.pi*5.3*MHz #Decay rate of |2> to |3>
gamma41 =  2*sc.pi*10.1*10**-9 #Decay rate of |4> to |1> used a lifetime of 31.2 s

#Laser Linewidths
gammalg = 0#2*sc.pi*2*MHz #493 laser linewidth
gammalr = 0#2*sc.pi*2*MHz #650 laser linewidth
gammalQ = 2*sc.pi*0.3002 #1762 laser linewidth

#Operators between |n> and |m> 
sig11 = basis(4,0) * basis(4,0).dag()
sig12 = basis(4,0) * basis(4,1).dag()
sig13 = basis(4,0) * basis(4,2).dag()
sig14 = basis(4,0) * basis(4,3).dag()
sig21 = basis(4,1) * basis(4,0).dag()
sig22 = basis(4,1) * basis(4,1).dag()
sig23 = basis(4,1) * basis(4,2).dag()
sig24 = basis(4,1) * basis(4,3).dag()
sig32 = basis(4,2) * basis(4,1).dag()
sig33 = basis(4,2) * basis(4,2).dag()
sig34 = basis(4,2) * basis(4,3).dag()
sig41 = basis(4,3) * basis(4,0).dag()
sig42 = basis(4,3) * basis(4,1).dag()
sig43 = basis(4,3) * basis(4,2).dag()
sig44 = basis(4,3) * basis(4,3).dag()

##Initial state of system  
psi0 = basis(4,0)

A = 2*sc.pi*3  # sweep rate
DeltaQinitial = -2*sc.pi*1 #Initial 1762 laser detuning
tstart = DeltaQinitial/A
tlist = np.linspace(tstart, -50.0*tstart, 100) #List of points for plotting purposes

##Decays and dissipations
C21 = np.sqrt(gamma21) * sig12 #Decay from |2> to |1>
C23 = np.sqrt(gamma23) * sig32 #Decay from |2> to |3>
C41 = np.sqrt(gamma41) * sig14 #Decay from |4> to |1>
Clg = np.sqrt(2*gammalg) * sig11 #From 493 laser linewidth
Clr = np.sqrt(2*gammalr) * sig33 #From 650 laser linewidth
ClQ = np.sqrt(2*gammalQ) * sig44 #From 1762 laser linewidth

singlePlot = 0

if singlePlot == 1:
    start_time = time.time()
    p_ex = qubit_integrate(Om14, DeltaQ, A, C41,ClQ, psi0, tlist)
    print('time elapsed = ' + str(time.time() - start_time))
    LZ= 1 - np.exp((-2*(np.pi)**2 * (Om14 **2)) / (A))
    print('Landau-Zener Aprroximation: ' + str(LZ))
    print('Hamiltonian Evolution: '+ str(np.real(p_ex)[-1]))
    
    fig, ax = plt.subplots(figsize=(12,8))
    ax.plot(tlist, np.real(p_ex), 'b', tlist, np.real(1-p_ex), 'r')
    ax.plot(tlist, LZ * np.ones(shape(tlist)), 'g')
    ax.set_xlabel('Time')
    ax.set_ylabel('Occupation probability')
    ax.set_title('Landau-Zener transition')
    ax.legend(("Excited state", "Ground state", "Landau-Zener formula"), loc=0);

if singlePlot == 0:
    yr = []
    xr = []
    A = 2*sc.pi*1  # sweep rate
    start_time = time.time()
    while A < 2*sc.pi*1000.3:
        tstart = DeltaQinitial/A
        tlist = np.linspace(tstart, -50*tstart, 100) #List of points for plotting purposes
        p_ex = qubit_integrate(Om14, DeltaQ, A, C41,ClQ, psi0, tlist)
        yr.append(np.real(p_ex)[-1])
        xr.append(A/(2*sc.pi))
        A += 2*sc.pi*0.1*A
    print('time elapsed = ' + str(time.time() - start_time))    
    fig, ax = plt.subplots(figsize=(12,8))
    ax.semilogx(xr, yr, 'ro',)
    ax.set_xlabel('Sweep Rate')
    ax.set_ylabel('Transfer Efficiency')
    ax.set_title('Landau-Zener transition');