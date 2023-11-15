# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 10:39:47 2023

@author: makdi
"""

import matplotlib.pyplot as plt

import scipy.constants as sc
import numpy as np
from qutip import *
from qutip import basis, steadystate, expect
qutip.settings.has_mkl=False

# initialise plot and line
# Ignore this
yr = []
xr = []
yg = []
xg = []
x = []
y = []
plt.ion()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, y, 'b-') 




#Start Here

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

#
#
#
#INPUT VALUES YOU WANT IN THIS SECTION
#
#
#

#Input Initial Values Here:
    #Units of frequnecies are in MZ
    #Time units will be in us
#Detunings
DeltaBlue = (2*sc.pi)*5 #detuning of 493 laser
DeltaRed = (2*sc.pi)*5 #detuning of 650 laser
DeltaLong= (2*sc.pi)*5 #detuning of 1762 laser

#493 beams
#Omg = 2*sc.pi*5 #Rabi frequrency of 2S1/2 to 2P1/2 
BluePiLight = 2*sc.pi*5 #Rabi frequrency of |1> to |3> and |2> to |4>
BlueSigmaPlusLight = 2*sc.pi*5 #Rabi frequrency of |1> to |4>
BlueSigmaMinusLight = 2*sc.pi*5 #Rabi frequrency of |2> to |3>

#650 beams
#Omr = 2*sc.pi*5 #Rabi frequrency of 2P1/2 to 2D3/2
RedPiLight = 2*sc.pi*5 #Rabi frequrency of |6> to |3> and |7> to |4>
RedSigmaPlusLight = 2*sc.pi*5 #Rabi frequrency of |5> to |3> and |6> to |4>
RedSigmaMinusLight = 2*sc.pi*5 #Rabi frequrency of |7> to |3> and |8> to |4>

#1762 beams
#OmL = 2*sc.pi*5 #Rabi frequrency of 2S1/2 to 2D3/2
LongPiLight = 2*sc.pi*5 #Rabi frequrency of |6> to |1> and |7> to |2>
LongSigmaPlusLight = 2*sc.pi*5 #Rabi frequrency of |7> to |1> and |8> to |2>
LongSigmaMinusLight = 2*sc.pi*5 #Rabi frequrency of |6> to |2> and |5> to |1>
#Decay Rates
#Decay Rates
SBlue=0.25
SRed=1
SLong=1
GammaBlue = BluePiLight/SBlue #Decay rate of 2P1/2 to 2S1/2
GammaRed =  RedPiLight/SRed #Decay rate of 2P1/2 to 2D3/2
GammaLong= LongPiLight/SLong #Decay rate of 2S1/2 to 2D3/2
GammaLaserBlue = 2*sc.pi*0 #493 laser linewidth
GammaLaserRed = 2*sc.pi*0 #650 laser linewidth
GammaLaserLong= 2*sc.pi*0 #1762 laser linewidth
B = 60/10000 #B-field in Tesla
wB = 2*sc.pi*5 #((sc.value('Bohr magneton')*B)/(sc.hbar))/1000000 #Larmor frequency in 2pi*MHz
Alpha=30 #Degrees angle of polarization vector. (Sin(Alpha),0,Cos(Alpha))




#Operators between |n> and |m> sig(row-col)
sig11 = basis(8,0) * basis(8,0).dag() #Connects |1> with |1>
sig13 = basis(8,0) * basis(8,2).dag() #Connects |1> with |3>
sig14 = basis(8,0) * basis(8,3).dag() #Connects |1> with |4>
sig22 = basis(8,1) * basis(8,1).dag() #Connects |2> with |2>
sig23 = basis(8,1) * basis(8,2).dag() #Connects |2> with |3>
sig24 = basis(8,1) * basis(8,3).dag() #Connects |2> with |4>
sig31 = basis(8,2) * basis(8,0).dag() #Connects |3> with |1>
sig32 = basis(8,2) * basis(8,1).dag() #Connects |3> with |2>
sig33 = basis(8,2) * basis(8,2).dag() #Connects |3> with |3>
sig35 = basis(8,2) * basis(8,4).dag() #Connects |3> with |5>
sig36 = basis(8,2) * basis(8,5).dag() #Connects |3> with |6>
sig37 = basis(8,2) * basis(8,6).dag() #Connects |3> with |7>
sig38 = basis(8,2) * basis(8,7).dag() #Connects |3> with |8>
sig41 = basis(8,3) * basis(8,0).dag() #Connects |4> with |1>
sig42 = basis(8,3) * basis(8,1).dag() #Connects |4> with |2>
sig44 = basis(8,3) * basis(8,3).dag() #Connects |4> with |4>
sig45 = basis(8,3) * basis(8,4).dag() #Connects |4> with |5>
sig46 = basis(8,3) * basis(8,5).dag() #Connects |4> with |6>
sig47 = basis(8,3) * basis(8,6).dag() #Connects |4> with |7>
sig48 = basis(8,3) * basis(8,7).dag() #Connects |4> with |8>
sig53 = basis(8,4) * basis(8,2).dag() #Connects |5> with |3>
sig54 = basis(8,4) * basis(8,3).dag() #Connects |5> with |4>
sig55 = basis(8,4) * basis(8,4).dag() #Connects |5> with |5>
sig63 = basis(8,5) * basis(8,2).dag() #Connects |6> with |3>
sig64 = basis(8,5) * basis(8,3).dag() #Connects |6> with |4>
sig66 = basis(8,5) * basis(8,5).dag() #Connects |6> with |6>
sig73 = basis(8,6) * basis(8,2).dag() #Connects |7> with |3>
sig74 = basis(8,6) * basis(8,3).dag() #Connects |7> with |4>
sig77 = basis(8,6) * basis(8,6).dag() #Connects |7> with |7>
sig83 = basis(8,7) * basis(8,2).dag() #Connects |8> with |3>
sig84 = basis(8,7) * basis(8,3).dag() #Connects |8> with |4>
sig88 = basis(8,7) * basis(8,7).dag() #Connects |8> with |8>

#constants to be used in Hamiltonian from Clebsch-Gordan
R3 = (np.sqrt(3)); R2 = (np.sqrt(2)) ;R6 = (np.sqrt(6)); Sin=np.sin(np.pi*Alpha/180); Cos=np.cos(np.pi*Alpha/180)
if np.abs(Sin)<0.00001: Sin=0
if np.abs(Cos)<0.00001: Cos=0

#Individual Hamiltionian parts
#Diagonals
H11 = (DeltaBlue-wB)*sig11
H22 = (DeltaBlue+wB)*sig22
H33 = (-wB/3)*sig33
H44 = (wB/3)*sig44
H55 = (DeltaRed-(6*wB/5))*sig55
H66 = (DeltaRed-(2*wB/5))*sig66
H77 = (DeltaRed+(2*wB/5))*sig77
H88 = (DeltaRed+(6*wB/5))*sig88
H_Dag=H11+H22+H33+H44+H55+H66+H77+H88
#S to P parts
H13 = (BluePiLight/R3)*Cos*sig13 
H31 = (BluePiLight/R3)*Cos*sig31
H14 = (-BlueSigmaPlusLight/R3)*Cos*sig14
H41 = (-BlueSigmaPlusLight/R3)*Cos*sig41
H23 = ((-BlueSigmaMinusLight)/R3)*Sin*sig23
H32 = ((-BlueSigmaMinusLight)/R3)*Sin*sig32
H24 = (-BluePiLight/R3)*Cos*sig24
H42 = (-BluePiLight/R3)*Cos*sig42
H_PS=H13+H31+H14+H41+H23+H32+H24+H42
#P to D parts
H53 = (-RedSigmaPlusLight/2)*Sin*sig53 
H35 = (-RedSigmaPlusLight/2)*Sin*sig35
H63 = (-RedPiLight/R3)*Cos*sig63
H36 = (-RedPiLight/R3)*Cos*sig36
H64 = (-RedSigmaPlusLight/(2*R3))*Sin*sig64
H46 = (-RedSigmaPlusLight/(2*R3))*Sin*sig46
H73 = (RedSigmaMinusLight/(2*R3))*Sin*sig73
H37 = (RedSigmaMinusLight/(2*R3))*Sin*sig37
H74 = (-RedPiLight/R3)*Cos*sig74
H47 = (-RedPiLight/R3)*Cos*sig47
H84 = (RedSigmaMinusLight/2)*Sin*sig84
H48 = (RedSigmaMinusLight/2)*Sin*sig48
H_PD=H53+H35+H63+H36+H64+H46+H73+H37+H74+H47+H84+H48

#Hamiltonian of system RWA
H = H_Dag+H_PS+H_PD

Hf = ((DeltaBlue-wB)*sig11 + ((-2/R3)*BluePiLight)*sig13 + (DeltaBlue+wB)*sig22 + ((2/R3)*BluePiLight)*sig24 + ((-2/R3)*BluePiLight)*sig31 + (-wB/3)*sig33 + ((1j/R2)*RedSigmaPlusLight)*sig35 + ((2/R6)*RedPiLight)*sig36 + ((-1j/R6)*RedSigmaMinusLight)*sig37 + ((2/R3)*BluePiLight)*sig42 + (wB/3)*sig44 + ((1j/R6)*RedSigmaPlusLight)*sig46 + ((2/R6)*RedPiLight)*sig47 + ((-1j/R2)*RedSigmaMinusLight)*sig48 + ((-1j/R2)*RedSigmaPlusLight)*sig53 + (DeltaRed-(6*wB/5))*sig55 + ((2/R6)*RedPiLight)*sig63 + ((-1j/R6)*RedSigmaPlusLight)*sig64 + (DeltaRed-(2*wB/5))*sig66 + ((1j/R6)*RedSigmaMinusLight)*sig73 + ((2/R6)*RedPiLight)*sig74 + (DeltaRed+(2*wB/5))*sig77 + ((1j/R2)*RedSigmaMinusLight)*sig84 + (DeltaRed+(6*wB/5))*sig88)
print(H)
print('sep')
print(Hf)
print('diff')
print(Hf-H)



##Initial state of system
psi0 = basis(8,0)

##Decays and dissipations
C41 = np.sqrt((2/3)*GammaBlue) * sig14 #Decay from |4> to |1>
C42 = -np.sqrt((1/3)*GammaBlue) * sig24 #Decay from |4> to |2>
C32 = np.sqrt((2/3)*GammaBlue) * sig23 #Decay from |3> to |2>
C31 = np.sqrt((1/3)*GammaBlue) * (sig13) #Decay from |3> to |1>
C35 = (np.sqrt((1/2)*GammaRed) * sig53) #Decay from |3> to |5> 
C36 = np.sqrt((1/3)*GammaRed) * (sig63) #Decay from |3> to |6>
C37 = (np.sqrt((1/6)*GammaRed) * sig73) #Decay from |3> to |7> 
C46 = (np.sqrt((1/6)*GammaRed) * sig64)# Decay from |4> to |6> 
C47 = np.sqrt((1/3)*GammaRed) * sig74 #Decay from |4> to |7>
C48 = (np.sqrt((1/2)*GammaRed) * sig84) #Decay from |4> to |8> 
Clg = np.sqrt(2*GammaLaserBlue) * (sig11 + sig22) #From 493 laser linewidth
Clr = np.sqrt(2*GammaLaserRed) * (sig55 + sig66 + sig77 + sig88) #From 650 laser linewidth
c_ops = [C41,C42,C32,C31,C35,C36,C37,C46,C47,C48,Clg,Clr]


DeltaRed =  -2*sc.pi*60
while DeltaRed < 2*sc.pi*40:
    DeltaBlue =  -2*sc.pi*10
    #H = ((DeltaBlue-wB)*sig11 + ((-2/R3)*BluePiLight)*sig13 + (DeltaBlue+wB)*sig22 + ((2/R3)*BluePiLight)*sig24 + ((-2/R3)*BluePiLight)*sig31 + (-wB/3)*sig33 + ((1j/R2)*RedSigmaPlusLight)*sig35 + ((2/R6)*RedPiLight)*sig36 + ((-1j/R6)*RedSigmaMinusLight)*sig37 + ((2/R3)*BluePiLight)*sig42 + (wB/3)*sig44 + ((1j/R6)*RedSigmaPlusLight)*sig46 + ((2/R6)*RedPiLight)*sig47 + ((-1j/R2)*RedSigmaMinusLight)*sig48 + ((-1j/R2)*RedSigmaPlusLight)*sig53 + (DeltaRed-(6*wB/5))*sig55 + ((2/R6)*RedPiLight)*sig63 + ((-1j/R6)*RedSigmaPlusLight)*sig64 + (DeltaRed-(2*wB/5))*sig66 + ((1j/R6)*RedSigmaMinusLight)*sig73 + ((2/R6)*RedPiLight)*sig74 + (DeltaRed+(2*wB/5))*sig77 + ((1j/R2)*RedSigmaMinusLight)*sig84 + (DeltaRed+(6*wB/5))*sig88)
    HTest = (DeltaBlue-wB)*sig11+(DeltaBlue+wB)*sig22+(-wB/3)*sig33+(wB/3)*sig44+(DeltaRed-(6*wB/5))*sig55+(DeltaRed-(2*wB/5))*sig66+(DeltaRed+(2*wB/5))*sig77+(DeltaRed+(6*wB/5))*sig88+(BluePiLight/R3)*Cos*sig13+(BluePiLight/R3)*Cos*sig31+(-BlueSigmaPlusLight/R3)*Cos*sig14+(-BlueSigmaPlusLight/R3)*Cos*sig41+((-BlueSigmaMinusLight)/R3)*Sin*sig23+((-BlueSigmaMinusLight)/R3)*Sin*sig32+(-BluePiLight/R3)*Cos*sig24+(-BluePiLight/R3)*Cos*sig42+(-RedSigmaPlusLight/2)*Sin*sig53+(-RedSigmaPlusLight/2)*Sin*sig35+(-RedPiLight/R3)*Cos*sig63+(-RedPiLight/R3)*Cos*sig36+(-RedSigmaPlusLight/(2*R3))*Sin*sig64+(-RedSigmaPlusLight/(2*R3))*Sin*sig46+(RedSigmaMinusLight/(2*R3))*Sin*sig73+(RedSigmaMinusLight/(2*R3))*Sin*sig37+(-RedPiLight/R3)*Cos*sig74+(-RedPiLight/R3)*Cos*sig47+(RedSigmaMinusLight/2)*Sin*sig84+(RedSigmaMinusLight/2)*Sin*sig48
    final_state = steadystate(HTest, c_ops) #Solve Hamiltonian for t = infinity
    fexpt1 = expect(sig44, final_state)  #Gives expectation value of solved Hamiltonian for excited state
    fexpt2 = expect(sig33, final_state)
    FexptPtot = 100*(fexpt1 + fexpt2)
    #print(FexptPtot)
    yr.append(FexptPtot)
    xr.append(DeltaRed/(2*sc.pi))
    DeltaRed += 2*sc.pi*0.3
    
plt.plot(xr, yr, 'r-')
plt.xlabel('ΔRed [MHz]')
plt.ylabel('Population in |3>+|4>')
plt.title('Alpha='+ str(Alpha)+', SBlue='+str(SBlue)+', SRed='+str(SRed)+', and Blue Detuning='+str(DeltaBlue/(2*sc.pi))+'MHz')
plt.show()