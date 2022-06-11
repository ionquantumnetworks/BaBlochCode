# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 18:48:00 2018

@author: James
"""
import matplotlib.pyplot as plt
import h5py

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
Deltag = -2*sc.pi*20 #detuning of 493 laser
Deltar = 2*sc.pi*10#detuning of 650 laser
DeltarEx = -2*sc.pi*(5)#detuning of 650 excite laser
#493 beams
Omgpi = 2*sc.pi*15#/(sqrt(100000)) #Rabi frequrency of |1> to |3> and |2> to |4>
Omgsp = 2*sc.pi*0#Rabi frequrency of |1> to |4>
Omgsm = 2*sc.pi*0 #Rabi frequrency of |2> to |3> 
#650 beams
Omrpi = 2*sc.pi*7#/(sqrt(100000))#Rabi frequrency of |6> to |3> and |7> to |4>
Omrsp = 2*sc.pi*7#/(sqrt(100000)) #Rabi frequrency of |5> to |3> and |6> to |4>
Omrsm = 2*sc.pi*15#(np.sqrt(1000))#Rabi frequrency of |7> to |3> and |8> to |4> 
#detunings
gammag =  2*sc.pi*15.1 #Decay rate of 2P1/2 to 2S1/2
gammar =  2*sc.pi*5.3 #Decay rate of 2P1/2 to 2D3/2
gammalg = 2*sc.pi*3#493 laser linewidth
gammalr = 2*sc.pi*0#650 laser linewidth
B = 5.23/10000 #B-field in Tesla
wB = ((sc.value('Bohr magneton')*B)/(sc.hbar))/1000000 #Larmor frequency in 2pi*MHz Bohr mag = 9.274*10^-24 J/T
#polarisation mising
#in terms of Rabi (ie 0.1: the bad polarisation results in the wrong transition to occur with 0.1 Rabi of the correct one)
mixc = 0.00 #mixing in the cleanupo pulse
mixp = 0.05 #mixing in the prep pulse
mixe = 0.0#mixing in the excite up pulse

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

#operators to input into mesolve - will yeild expectation value of these operators
e_ops = [sig11,sig22] #S-levels
#e_ops = [sig33,sig44] #P-levels
#e_ops = [sig55,sig66,sig77,sig88] #D-levels

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
#H0 = ((Deltag-wB)*sig11 + ((-2/R3)*0)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*0)*sig24 + ((-2/R3)*0)*sig31 + (-wB/3)*sig33 + ((1j/R2)*0)*sig35 + ((2/R6)*0)*sig36 + ((-1j/R6)*0)*sig37 + ((2/R3)*0)*sig42 + (wB/3)*sig44 + ((1j/R6)*0)*sig46 + ((2/R6)*0)*sig47 + ((-1j/R2)*0)*sig48 + ((-1j/R2)*0)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*0)*sig63 + ((-1j/R6)*0)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*0)*sig73 + ((2/R6)*0)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*0)*sig84 + (DeltarEx+(6*wB/5))*sig88)
Hp = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*mixp*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*mixp*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*mixp*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*mixp*Omrsm)*sig84 + (DeltarEx+(6*wB/5))*sig88)
Hp2 = ((Deltag-wB)*sig11 + ((-2/R3)*0)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*0)*sig24 + ((-2/R3)*0)*sig31 + (-wB/3)*sig33 + ((1j/R2)*0)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*(mixc*Omrsm))*sig37 + ((2/R3)*0)*sig42 + (wB/3)*sig44 + ((1j/R6)*0)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*(mixc*Omrsm))*sig48 + ((-1j/R2)*0)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*0)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)**(mixc*Omrsm))*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*(mixc*Omrsm)*sig84) + (DeltarEx+(6*wB/5))*sig88)
He = ((Deltag-wB)*sig11 + ((-2/R3)*0)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*0)*sig24 + ((-2/R3)*0)*sig31 + (-wB/3)*sig33 + ((1j/R2)**mixe*Omrsm)*sig35 + ((2/R6)*0)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*0)*sig42 + (wB/3)*sig44 + ((1j/R6)*mixe*Omrsm)*sig46 + ((2/R6)*0)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*mixe*Omrsm)*sig53 + (DeltarEx-(6*wB/5))*sig55 + ((2/R6)*0)*sig63 + ((-1j/R6)*mixe*Omrsm)*sig64 + (DeltarEx-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*0)*sig74 + (DeltarEx+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (DeltarEx+(6*wB/5))*sig88)
#He = ((Deltag-wB)*sig11 + ((-2/R3)*0)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*0)*sig24 + ((-2/R3)*0)*sig31 + (-wB/3)*sig33 + ((1j/R2)*0)*sig35 + ((2/R6)*0)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*0)*sig42 + (wB/3)*sig44 + ((1j/R6)*0)*sig46 + ((2/R6)*0)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*0)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*0)*sig63 + ((-1j/R6)*0)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*0)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (DeltarEx+(6*wB/5))*sig88)
#Ht = ((Deltag-wB)*sig11 + ((-2/R3)*Omgpi)*sig13 + (Deltag+wB)*sig22 + ((2/R3)*Omgpi)*sig24 + ((-2/R3)*Omgpi)*sig31 + (-wB/3)*sig33 + ((1j/R2)*Omrsp)*sig35 + ((2/R6)*Omrpi)*sig36 + ((-1j/R6)*Omrsm)*sig37 + ((2/R3)*Omgpi)*sig42 + (wB/3)*sig44 + ((1j/R6)*Omrsp)*sig46 + ((2/R6)*Omrpi)*sig47 + ((-1j/R2)*Omrsm)*sig48 + ((-1j/R2)*Omrsp)*sig53 + (Deltar-(6*wB/5))*sig55 + ((2/R6)*Omrpi)*sig63 + ((-1j/R6)*Omrsp)*sig64 + (Deltar-(2*wB/5))*sig66 + ((1j/R6)*Omrsm)*sig73 + ((2/R6)*Omrpi)*sig74 + (Deltar+(2*wB/5))*sig77 + ((1j/R2)*Omrsm)*sig84 + (DeltarEx+(6*wB/5))*sig88)
#H = H11+H22+H33+H44+H55+H66+H77+H88+H13+H31+H14+H41+H23+H32+H24+H42+H53+H35+H54+H45+H64+H46+H73+H37+H74+H47+H83+H38+H84+H48

#Time dependence of the drive
def Ht_coeffPrep(t,args):
    on = -15.0
    pulse = 8
    edge = 150
    off = on + pulse
    y4 = 0.5*(np.tanh(edge*(t-on))-np.tanh(edge*(t-off)))
    return y4

def Ht_coeffPrep2(t,args):
    on = -6
    pulse = 1
    edge = 150
    off = on + pulse
    y4 = 0.5*(np.tanh(edge*(t-on))-np.tanh(edge*(t-off)))
    return y4

def Ht_coeffEx(t,args):
    on = 0.994
    pulse = 0.2
    edge = 200
    off = on + pulse
    y4 = 0.5*(np.tanh(edge*(t-on))-np.tanh(edge*(t-off)))
    return y4


#Complete Hamiltonian
H = [[Hp, Ht_coeffPrep],[Hp2,Ht_coeffPrep2],[He, Ht_coeffEx]] #[Hp2,Ht_coeffPrep2]

#Initial state of system
#psi0 = (1/(np.sqrt(4)))*(basis(8,4)+basis(8,5)+basis(8,6)+basis(8,7)) #D-state superposition
#psi0 = (1/(np.sqrt(8)))*(basis(8,0)+basis(8,1)+basis(8,2)+basis(8,3)+basis(8,4)+basis(8,5)+basis(8,6)+basis(8,7))
psi0 = (1/(np.sqrt(2)))*(basis(8,0)+basis(8,1)) #ground state superposition
#psi0 = basis(8,7) #stretch-state
#psi0 = (1/(np.sqrt(2)))*(basis(8,2)+basis(8,3)) #P-state superposition
#psi0 = basis(8,4)
#psi0 = np.sqrt(0.05)*basis(8,0)+np.sqrt(0.05)*basis(8,1)+np.sqrt(0.03)*basis(8,2)+np.sqrt(0.03)*basis(8,3)+np.sqrt(0.22)*basis(8,4)+np.sqrt(0.17)*basis(8,5)+np.sqrt(0.20)*basis(8,6)+np.sqrt(0.25)*basis(8,7)
#psi0 = (1/(np.sqrt(3)))*(basis(8,5)+basis(8,6)+basis(8,7)) #D-state superposition
#STEADY STATE POPS: 0.04797718762193412 0.04662544116872065 0.03470595961939421 0.034370105505866694 0.2170933242442322 0.17455678227670463 0.19925282499341132 0.24541837456973623
#psi0 = (np.sqrt(0.015)*basis(8,4)+ np.sqrt(0.014)*basis(8,5)+ np.sqrt(0.013)*basis(8,6)+ np.sqrt(0.86)*basis(8,7)) #D-state superposition

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


startT = -15
stopT = 1.1 #stop time for calculating fidelity
points = (stopT - startT)*1000 #ensures 1 ns point spaceing
times = np.linspace(startT,stopT, int(points))
result = mesolve(H, psi0, times, c_ops, [sig11,sig22,sig33,sig44,sig55,sig66,sig77,sig88])
#summs all the data over all time
allCounts = result.expect[2]+result.expect[3]
allGood = result.expect[3]
allBad = result.expect[2]
windowStart = int(points)-500 #should be at 900 ns
windowWidth = 500
windowStop = windowStart + windowWidth
#pulls out and sums the last 100 points - as each point is set up to be 1 ns this is 900 ns - 1300 ns corrisponding to the photon
totalPhoton = sum(allCounts[windowStart:windowStop])
goodPhoton = sum(allGood[windowStart:windowStop])
badPhoton =sum(allBad[windowStart:windowStop])
#scale factor to scale qutip data to real-life data
photpwind = (gammag*10**-3) #expected 493 photons per ns given a population = 1 (*10**6 for gamma in ns is *10**-3)
br = 0.75 #branching
col = 0.08 #0.6 NA into fiber
qe = 0.4 #APD
coup = 0.20 #estimated coupling
runs = 557634 #total number of experiment runs
scale = (runs/1190)#photpwind*totalPhoton*runs*coup*qe*col*0.75# (br*col*qe*coup*runs)/totalPhoton

#various plotting
fig, ax = subplots()
ax.plot((result.times)*1000, (result.expect[0]+result.expect[1]));#Ground State
ax.plot((result.times)*1000, (result.expect[2]+result.expect[3]));#P-levels
ax.plot((result.times)*1000, (result.expect[4]+result.expect[5]+result.expect[6]+result.expect[7]));#D-levels
ax.set_xlabel('Time [ns]');
ax.set_ylabel('Population');
ax.legend(("S","P","D"));
show()
fig, ax = subplots()
ax.plot((result.times)*1000, (result.expect[4]));
ax.plot((result.times)*1000, (result.expect[5]));
ax.plot((result.times)*1000, (result.expect[6]));
ax.plot((result.times)*1000, (result.expect[7]));
ax.set_xlabel('Time [ns]');
ax.set_ylabel('Population');
ax.legend(("5","6","7","8"));
show()
fig, ax = subplots()
ax.plot((result.times)*1000, (result.expect[2]));
ax.plot((result.times)*1000, (result.expect[3]));
ax.set_xlabel('Time [ns]');
ax.set_ylabel('Population');
ax.legend(("3","4"));
show()
Photon = ((result.expect[2]+result.expect[3]))*0.75
fig, ax = subplots()
ax.plot((result.times)*1000,Photon);#P-levels
ax.set_xlabel('Time [ns]');
ax.set_ylabel('');
ax.legend(("Photon Shape",""));
show()

print("S:" + str(result.expect[0][-1]+result.expect[1][-1]))
print("1:" + str(result.expect[0][-1]))
print("2:" + str(result.expect[1][-1]))
print("P:" + str(result.expect[2][-1]+result.expect[3][-1]))
print("D:" + str(result.expect[4][-1]+result.expect[5][-1]+result.expect[6][-1]+result.expect[7][-1]))
print("5:" + str(result.expect[4][-1]))
print("6:" + str(result.expect[5][-1]))
print("7:" + str(result.expect[6][-1]))
print("8:" + str(result.expect[7][-1]))

#Save graph data
#output_data = np.vstack((times*1000, (result.expect[2]+result.expect[3]))) # join time and expt˓→data
#file_data_store('E:\\IonTrapData\\DPrep Photon Shapes\\493PhotonShape.dat', output_data.T, numtype="real") # Note the .T for transpose!
#thefile = open('E:\\IonTrapData\\DPrep Photon Shapes\\493PhotonShape.dat', 'w')
#for item in Photon:
#  thefile.write("%s\n" % item)
#thefile.close()

with h5py.File('000007129-PhotonShapeTempMeasure.h5', 'r') as f:
    # List all groups that could be imported
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[0]
    #read data
    ClickDeltas = list(f['datasets/ClickDeltas'])
    info = np.array(f.get('expid'))
    print(info)

#calculate max counts and calculate how many bins the plot will have
maxn = round(max(ClickDeltas))
minn = 0#round(min(ClickDeltas))
binsg = round(3000)
#histogram data
TimeData, y1 = np.histogram(ClickDeltas, bins=np.linspace(minn, maxn, binsg))
print(maxn-minn)
#plots histogram of data
fig, ax1 = plt.subplots()
ax1.stairs(TimeData,y1,hatch='//')
ax1.plot((result.times)*1000,Photon*scale)
ax1.plot((result.times)*1000, (0.75*result.expect[2])*scale)
ax1.plot((result.times)*1000, (0.75*result.expect[3])*scale)
ax1.plot(times*1000,Ht_coeffEx(times,1)*100)
ax1.set_xlim([990,1150])
ax1.set_ylim([0,120])
ax1.set_ylabel('Occurances')
ax1.set_xlabel('Time (ns)')

#photon analysis
FGP = goodPhoton/totalPhoton
FBP = badPhoton/totalPhoton
print('Fraction Good Photon ' + str(FGP))
print('Fraction Bad Photon ' + str(FBP))

#Fidelities
Target = (basis(4,0)+basis(4,3)).unit()
TargetDensity = Target*Target.dag()
#print(Target)

Psi1 = basis(4,0)
Psi2 = basis(4,1)
Psi3 = basis(4,2)
Psi4 = basis(4,3)

badState = (basis(4,1)+basis(4,2)).unit()
badDensity = badState*badState.dag()
FinalDensity = (FGP*TargetDensity)+(FBP*badDensity)

FidFinal = Target.dag()*(FinalDensity)*Target
print(FidFinal)

N =2
lbls_list = [[str(d) for d in range(N)], ["u", "d"]]

xlabels = []

for inds in tomography._index_permutations([len(lbls) for lbls in lbls_list]):
     xlabels.append("".join([lbls_list[k][inds[k]]
                            for k in range(len(lbls_list))]))

fig, ax = matrix_histogram(FinalDensity,xlabels, xlabels) # xlabels=xlabels, ylabels=xlabels)
plt.show()
