# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 18:48:00 2018

@author: James
"""
import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
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

#Units of frequnecies are in MZ
#Time units will be in us
Deltag = 0*sc.pi*0 #detuning of 493 laser
Deltar = 0*sc.pi*0 #detuning of 650 laser
#493 beams
#Omg = 2*sc.pi*5 #Rabi frequrency of 2S1/2 to 2P1/2 
Omgpi = 2*sc.pi*10 #Rabi frequrency of |1> to |3> and |2> to |4>
Omgsp = 2*sc.pi*10#Rabi frequrency of |1> to |4>
Omgsm = 2*sc.pi*10 #Rabi frequrency of |2> to |3>
#650 beams
#Omr = 2*sc.pi*5 #Rabi frequrency of 2P1/2 to 2D3/2
Omrpi = 2*sc.pi*5 #Rabi frequrency of |6> to |3> and |7> to |4>
Omrsp = 2*sc.pi*5 #Rabi frequrency of |5> to |3> and |6> to |4>
Omrsm = 2*sc.pi*5 #Rabi frequrency of |7> to |3> and |8> to |4>
#detunings
gammag =  2*sc.pi*15.1 #Decay rate of 2P1/2 to 2S1/2
gammar =  2*sc.pi*5.3 #Decay rate of 2P1/2 to 2D3/2
gammalg = 2*sc.pi*0 #493 laser linewidth
gammalr = 2*sc.pi*0 #650 laser linewidth
B = 5/10000 #B-field in Tesla
#tlist = np.linspace(0, 0.25, 10) #List of points for plotting purposes
wB = ((sc.value('Bohr magneton')*B)/(sc.hbar))/1000000 #Larmor frequency in 2pi*MHz

#Operators between |n> and |m> sig(row-col)
sig11 = basis(8,0) * basis(8,0).dag()
sig13 = basis(8,0) * basis(8,2).dag()
sig14 = basis(8,0) * basis(8,3).dag()
sig22 = basis(8,1) * basis(8,1).dag()
sig23 = basis(8,1) * basis(8,2).dag()
sig24 = basis(8,1) * basis(8,3).dag()
sig31 = basis(8,2) * basis(8,0).dag()
sig32 = basis(8,2) * basis(8,1).dag()
sig33 = basis(8,2) * basis(8,2).dag()
sig35 = basis(8,2) * basis(8,4).dag()
sig36 = basis(8,2) * basis(8,5).dag()
sig37 = basis(8,2) * basis(8,6).dag()
sig38 = basis(8,2) * basis(8,7).dag()
sig41 = basis(8,3) * basis(8,0).dag()
sig42 = basis(8,3) * basis(8,1).dag()
sig44 = basis(8,3) * basis(8,3).dag()
sig45 = basis(8,3) * basis(8,4).dag()
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

e_ops1 = [sig44] #operators to input into mesolve
e_ops2 = [sig33]

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
H31 = (Omgpi/R3)*sig31
H14 = (-Omgsp/R3)*sig14
H41 = ((1j*Omgsp)/R3)*sig41
H23 = ((-Omgsm)/R3)*sig23
H32 = ((-1j*Omgsm)/R3)*sig32
H24 = (-Omgpi/R3)*sig24
H42 = (-Omgpi/R3)*sig42
#P to D parts
H53 = (-Omrsp/2)*sig53 
H35 = ((Omrsp*1j)/2)*sig35 
H54 = 0*sig54
H45 = 0*sig45
H63 = (-Omrpi/R3)*sig63
H36 = (-Omrpi/R3)*sig36
H64 = (-Omrsp/(2*R3))*sig64
H46 = ((Omrsp*1j)/(2*R3))*sig46
H73 = (Omrsm/(2*R3))*sig73
H37 = ((Omrsm*1j)/(2*R3))*sig37
H74 = (-Omrpi/R3)*sig74
H47 = (-Omrpi/R3)*sig47
H83 = 0*sig83
H38 = 0*sig38
H84 = (Omrsm/2)*sig84
H48 = ((Omrsm*1j)/2)*sig48
#Hamiltonian of system RWA
H = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (Deltar+(6*wB/5))*sig88)
#H = H11+H22+H33+H44+H55+H66+H77+H88+H13+H31+H14+H41+H23+H32+H24+H42+H53+H35+H54+H45+H64+H46+H73+H37+H74+H47+H83+H38+H84+H48

##Initial state of system
psi0 = basis(8,0)

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


#Solutions
#x = 1000000000
#n = mesolve(H, psi0, tlist, c_ops, e_ops, options=Options(nsteps=x)) #Solves Hamiltionian at point on tlist
#final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
#fexpt1 = expect(e_ops1, final_state) #Calculates expectation values (e_ops list) for Hamiltonion at t = inifinity
#fexpt2 = expect(e_ops2, final_state) #Calculates expectation values (e_ops list) for Hamiltonion at t = inifinity
#print(fexpt1,fexpt2,fexpt1+fexpt2)
#plot_expectation_values(n, show_legend=True,figsize=(8, 8))

yr = []
xr = []
yg = []
xg = []
x = []
y = []

# initialise plot and line
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'b-') 

Deltag = -2*sc.pi*100
while Deltag < 2*sc.pi*50:
    Deltar = -2*sc.pi*29 #detuning of 650 laser
    H = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (Deltar+(6*wB/5))*sig88)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt1 = expect(sig44, final_state)  #Gives expectation value of solved Hamiltonian for excited state
    fexpt2 = expect(sig33, final_state)
    FexptPtot = fexpt1 + fexpt2
    yg.append(FexptPtot)
    xg.append(Deltag/(2*sc.pi))
    Deltag += 2*sc.pi*1

plt.plot(xg, yg,'c-')
plt.xlabel('Δg [MHz]')
plt.ylabel('Population in |3>+|4>')
plt.show()

Deltar =  -2*sc.pi*60
while Deltar < 2*sc.pi*40:
    Deltag =  -2*sc.pi*99
    H = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (Deltar+(6*wB/5))*sig88)
    final_state = steadystate(H, c_ops) #Solve Hamiltonian for t = infinity
    fexpt1 = expect(sig44, final_state)  #Gives expectation value of solved Hamiltonian for excited state
    fexpt2 = expect(sig33, final_state)
    FexptPtot = fexpt1 + fexpt2
    #print(FexptPtot)
    yr.append(FexptPtot)
    xr.append(Deltar/(2*sc.pi))
    Deltar += 2*sc.pi*1
    
plt.plot(xr, yr, 'r-')
plt.xlabel('Δr [MHz]')
plt.ylabel('Population in |3>+|4>')
plt.show()

print(wB)
##thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\493 Photon Shape Tests\\DeltaR_8Mar.dat', 'w')
##for item in xr:
##  thefile.write("%s\n" % item)
##thefile.close()

##thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\493 Photon Shape Tests\\PopPR_8Mar.dat', 'w')
##for item in yr:
##  thefile.write("%s\n" % item)
##thefile.close()
###
#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\493 Photon Shape Tests\\DeltaR493-8Mar.dat', 'w')
#for item in xg:
#  thefile.write("%s\n" % item)
#thefile.close()
##
#thefile = open('G:\\Team Drives\\Ions\\03 - Projects\\Current Projects\\Rb Ba+ hybrid\\493 Photon Shape Tests\\PopP493-8Mar.dat', 'w')
#for item in yg:
#  thefile.write("%s\n" % item)
#thefile.close()