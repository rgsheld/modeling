# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 19:37:07 2017

@author: Sergei
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 100  # number of release sites
PoccI = 1  # probability of a site being occupied at the start of the sim
N_occ = N * PoccI
PrI = 0.8  # probability of a occupied site releasing a vesicle at the start
Tau_f = 6   # time constant of facilitation in stim number
delta_P = 0.1  # release probability increment due to facilitation
K_f = 0.25  # rate of replenishment per available site / stim
Tau_R = 3  # temporal delay of replenishment
K_ud = 0.0001 # rate of vesicle undocking / stim
K_r = PrI + K_ud # total rate of site opening / stim


# Initialize data arrays
Release = np.zeros([50,1])
Kon = np.zeros([50,1])
Pr = np.zeros([50,1])

# 1st stimulus
Release[0,:] = (N_occ*PrI) 
Kon[0,:] = K_f*(N-N_occ)
Pr[0,:] = PrI

# Update all the parameters
Kf_plus = K_f*(N - N_occ)
K_f = (Kf_plus) + ((1 - np.exp(-1/Tau_R))*(Kon[0, :]-Kf_plus))
N_occ = N_occ + K_f*(N - N_occ)*PrI - K_r * N_occ
#N_occ = (N_occ - Release[0, :]) + K_f*(N -(N_occ - Release[0, :]))  # release site occupancy with scaling replenishment
PrN_plus = PrI + delta_P*(1-PrI) # release prob immediately after last AP
PrN = (PrN_plus) + ((1 - np.exp(-1/Tau_f))*(PrI-PrN_plus))

Release[1,:] = (N_occ*PrN)
Kon[1,:] = K_f*(N -(N_occ - Release[0, :]))
Pr[1,:] = PrN

for i in range(2, 50):      
    N_occ = N_occ + K_f*(N - N_occ)*PrN - K_r * N_occ
    #N_occ = (N_occ - Release[i, :]) + K_f*(N -(N_occ - Release[i, :]))
    PrN_plus = PrN + delta_P*(1-PrN)
    PrN = (PrN_plus) + ((1 - np.exp(-1/Tau_f))*(PrI-PrN_plus))

    Release[i,0] = (N_occ*PrN)
    Kon[i,:] = K_f*(N -(N_occ - Release[i, :]))
    Pr[i,:] = PrN


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True)
ax1.plot(Release, 'ro-')
ax1.set_title('Amplitude')
ax2.plot(Pr, 'ro-')
ax2.set_title('Release Probability')
ax3.plot(Kon, 'ro-')
ax3.set_title('Replenishment')