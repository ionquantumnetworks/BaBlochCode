# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 12:30:38 2019

@author: James
"""

import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from pylab import *
from qutip import *
qutip.settings.has_mkl=False

#       ----- mj=+1/2 |4>
# 2P1/2 
#       ----- mj=-1/2 |3>
#                               ----- mj=+3/2 |8>
#                               ----- mj=+1/2 |7>
#                         2D3/2
#                               ----- mj=-1/2 |6>
#                               ----- mj=-3/2 |5>
#       ----- mj=+1/2 |2>
# 2S1/2 
#       ----- mj=-1/2 |1>

#Units of frequnecies are in MHz
#Time units will be in us
#Detunigs
Deltag = -2*sc.pi*20 #detuning of 493 laser
Deltar = 2*sc.pi*0#detuning of 650 laser
#493 beams
Omgpi = 2*sc.pi*20#Rabi frequrency of |1> to |3> and |2> to |4>
Omgsp = 2*sc.pi*0 #Rabi frequrency of |1> to |4>
Omgsm = 2*sc.pi*0 #Rabi frequrency of |2> to |3>
#650 beams
Omrpi = 2*sc.pi*00 #Rabi frequrency of |6> to |3> and |7> to |4>
Omrsp = 2*sc.pi*5 #Rabi frequrency of |5> to |3> and |6> to |4>
Omrsm = 2*sc.pi*5 #Rabi frequrency of |7> to |3> and |8> to |4> 
#Linewidths
gammag =  2*sc.pi*15.1 #Decay rate of 2P1/2 to 2S1/2
gammar =  2*sc.pi*5.3 #Decay rate of 2P1/2 to 2D3/2
gammalg = 2*sc.pi*3 #493 laser linewidth
gammalr = 2*sc.pi*3 #650 laser linewidth
B = 5.23/10000 #B-field in Tesla
#tlist = np.linspace(0, 0.06, 2000) #List of points for plotting purposes
wB = ((sc.value('Bohr magneton')*B)/(sc.hbar))/1000000 #Larmor frequency in 2pi*MHz Bohr mag = 9.274*10^-24 J/T

#Operators between |n> and |m> sig(row-col)
sig11 = basis(8,0) * basis(8,0).dag()
sig13 = basis(8,0) * basis(8,2).dag()
sig14 = basis(8,0) * basis(8,3).dag()
sig22 = basis(8,1) * basis(8,1).dag()
sig23 = basis(8,1) * basis(8,2).dag()
sig24 = basis(8,1) * basis(8,3).dag()
sig31 = basis(8,2) * basis(8,0).dag()
sig33 = basis(8,2) * basis(8,2).dag()
sig35 = basis(8,2) * basis(8,4).dag()
sig36 = basis(8,2) * basis(8,5).dag()
sig37 = basis(8,2) * basis(8,6).dag()
sig42 = basis(8,3) * basis(8,1).dag()
sig44 = basis(8,3) * basis(8,3).dag()
sig46 = basis(8,3) * basis(8,5).dag()
sig47 = basis(8,3) * basis(8,6).dag()
sig48 = basis(8,3) * basis(8,7).dag()
sig53 = basis(8,4) * basis(8,2).dag()
sig54 = basis(8,4) * basis(8,3).dag()
sig55 = basis(8,4) * basis(8,4).dag()
sig63 = basis(8,5) * basis(8,2).dag()
sig64 = basis(8,5) * basis(8,3).dag()
sig66 = basis(8,5) * basis(8,5).dag()
sig73 = basis(8,6) * basis(8,2).dag()
sig74 = basis(8,6) * basis(8,3).dag()
sig77 = basis(8,6) * basis(8,6).dag()
sig83 = basis(8,7) * basis(8,2).dag()
sig84 = basis(8,7) * basis(8,3).dag()
sig88 = basis(8,7) * basis(8,7).dag()

#constants
R3 = (np.sqrt(3))
R2 = (np.sqrt(2))
R6 = (np.sqrt(6))

#operators to input into mesolve
#e_ops = [sig11,sig22] #S-levels
#e_ops = [sig33,sig44] #P-levels
#e_ops = [sig55,sig66,sig77,sig88] #D-levels
#e_ops33 = [sig33]
#e_ops44 = [sig44]

#Individual Hamiltionian parts
#Diagonals
H11 = (Deltag-wB)*sig11
H22 = (Deltag+wB)*sig22
H33 = (-wB/3)*sig33
H44 = (wB/3)*sig44
H55 = (Deltar-(6*wB/5))*sig55
H66 = (Deltar-(2*wB/5))*sig66
H77 = (Deltar+(2*wB/5))*sig77
H88 = (Deltar+(6*wB/5))*sig88
#S to P parts
H13 = (Omgpi/R3)*sig13
H14 = (-Omgsp/R3)*sig14
H23 = (-Omgsm/R3)*sig23
H24 = (-Omgpi/R3)*sig24
#P to D parts
H53 = (-Omrsp/2)*sig53 
H54 = 0*sig54
H63 = (-Omrpi/R3)*sig63
H64 = (-Omrsp/(2*R3))*sig64
H73 = (Omrsm/(2*R3))*sig73
H74 = (-Omrpi/R3)*sig74
H83 = 0*sig83
H84 = (Omrsm/2)*sig84
#Hamiltonian of system RWA
H = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (Deltar+(6*wB/5))*sig88)
#H = H11 + H22 + H33 + H44 + H55 + H66 + H77 + H88 + H13 + H14 + H23 + H24 + H53 + H54 + H64 + H73 + H74 + H83 + H84

#Initial state of system
#psi0 = (1/(np.sqrt(4)))*(basis(8,4)+basis(8,5)+basis(8,6)+basis(8,7)) #D-state superposition
psi0 = (1/(np.sqrt(2)))*(basis(8,0)+basis(8,1)) #ground state superposition
#psi0 = (1/(np.sqrt(2)))*(basis(8,2)+basis(8,3)) #excited state superposition
#psi0 = basis(8,7) #stretch-state

##Decays and dissipations
C41 = np.sqrt((2/3)*gammag) * sig14 #Decay from |4> to |1>
C42 = np.sqrt((1/3)*gammag) * sig24 #Decay from |4> to |2>
C32 = np.sqrt((2/3)*gammag) * sig23 #Decay from |3> to |2>
C31 = np.sqrt((1/3)*gammag) * (sig13) #Decay from |3> to |1>
C35 = (np.sqrt((1/2)*gammar) * sig53) #Decay from |3> to |5> 
C36 = np.sqrt((1/3)*gammar) * (sig63)
C37 = (np.sqrt((1/6)*gammar) * sig73) #Decay from |3> to |7> 
C46 = (np.sqrt((1/6)*gammar) * sig64)# Decay from |4> to |6> 
C47 = np.sqrt((1/3)*gammar) * sig74
C48 = (np.sqrt((1/2)*gammar) * sig84) #Decay from |4> to |8> 
Clg = np.sqrt(2*gammalg) * (sig11 + sig22) #From 493 laser linewidth
Clr = np.sqrt(2*gammalr) * (sig55 + sig66 + sig77 + sig88) #From 650 laser linewidth
c_ops = [C41,C42,C32,C31,C35,C36,C37,C46,C47,C48,Clg,Clr]

# initialise plot and line
#tlist = np.linspace(0,100,200) #List of points for plotting purposes  
#final_state = steadystate(H, c_ops)
#n = mesolve(H, final_state, tlist, c_ops,e_ops).expect[0]
#a_op = [basis(8,0)]
##tlist,g2=coherence_function_g2(H, final_state, tlist, c_ops, sig13)
#G2 = correlation_2op_1t(H, final_state, tlist, c_ops, sig13, sig11)
#g2 = (G2.real)/G2[-1].real
##(G2.real)//(n*n)G2[1].real
#plot(tlist,g2)
##plot(tlist,g2)
##show()
#print(G2[-1].real)


#plotting state evolution
times = np.linspace(0, 1.2, 2000)
result = mesolve(H, psi0, times, c_ops, [sig11,sig22,sig33,sig44,sig55,sig66,sig77,sig88])
fig, ax = subplots()
ax.plot((result.times)*1000, (result.expect[0]+result.expect[1]));#Ground State
ax.plot((result.times)*1000, (result.expect[2]+result.expect[3]));#P-levels
ax.plot((result.times)*1000, (result.expect[4]+result.expect[5]+result.expect[6]+result.expect[7]));#D-levels
ax.set_xlabel('Time [ns]');
ax.set_ylabel('Population');
ax.legend(("S","P","D"));
show()
fexpt3 = expect(sig33, final_state)
fexpt4 = expect(sig44, final_state)
PopinP = fexpt3 + fexpt4
print(PopinP)

#g2 by hand method
times = np.linspace(0, 0.030, 5000)
result = mesolve(H, psi0, times, c_ops, [sig11,sig22,sig33,sig44,sig55,sig66,sig77,sig88])
final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
fexpt4 = expect(sig44, final_state)  #Gives expectation value of solved Hamiltonian for excited state
fexpt3 = expect(sig33, final_state)
FexptPtot = fexpt4 + fexpt3
g2 = (result.expect[2]+result.expect[3])/FexptPtot


fig, ax = subplots()
ax.plot(-(result.times)*1000, g2,'b', (result.times)*1000, g2, 'b');
ax.set_xlabel('Time [ns]');
ax.set_ylabel('g2');
ax.legend(("g2", ""));
show()

#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Completed Projects\\780 1stage conversion g2\\Ba8lvlg2.dat', 'w')
#for item in g2:
#  thefile.write("%s\n" % item)
#thefile.close()


#print(H)
#n33 = mesolve(H, final_state, tlist, c_ops, e_ops33).expect[0]
#n44 = mesolve(H, final_state, tlist, c_ops, e_ops44).expect[0]
#final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
#fexpt1 = expect(sig44, final_state)  #Gives expectation value of solved Hamiltonian for excited state
#fexpt2 = expect(sig33, final_state)
#FexptPtot = fexpt1 + fexpt2
##plot_expectation_values(n33, show_legend=True,figsize=(8, 8))
#n3344 = (n33+n44)
##plot(tlist,n3344)
#g2 = n3344/(FexptPtot)
#plot(tlist,g2)
#print(g2)