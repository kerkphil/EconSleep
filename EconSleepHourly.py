# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 12:12:25 2017

@author: Kerk Phillips
"""

import numpy as np
import matplotlib.pyplot as plt

def moddefs(B, Z, y, A, C, *mparams):
    '''
    B is the stock of savings now
    Z is the stock of sleep capital now
    y is the value of the circadian cycle now
    A is an indicator for working, playing or sleeping now
    P is a binary indicating playing now
    C is consumption now if P=1
    
    UC is utility from consumption >= 0
    Uy is utility from following the circadian cycle
    Bp is savings next period
    Zp is sleep stock next period
    U is total utility
    '''
    # effective labor as a function of sleep stock
    # e = h*(Z**(1. - eta)+)/(1. - eta)
    e = h*Z**eta
    # upward shift in circadian cycle as function of sleep stock
    d = np.arctan(f*Z)/np.pi + g
    # if working (A=0)
    if (A == 0):
        UC = 0.
        Uy = y*d
        Bp = B + e*wbar
        Zp = Z*(1-delW)
    # if playing (A=1)
    elif (A == 1):
        # consumption greater than savings is not allowed
        if C > B:
            UC = -1.E+99
            Bp = 0.
        else:
            UC = e*C**gamma
            Bp = B - C
        Uy = y*d
        Zp = Z*(1-delP)
    else:
        UC = 0.
        Uy = -y*d
        Bp = B
        Zp = Z + 1/S
    # total utility
    U = UC + chiS*Uy
    return U, Bp, Zp


def nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

# main program

# define parameters
S = 24        # number of periods per day
eta = .5      # curvature of phi(Z) functuion
h = 1.        # scaling factor for phi(Z) function
wbar = 1.     # constant wage
f = .5        # scaling factor for d(Z) function
g = .1        # additive factor for d(Z) function
b = 4.        # phase shifter for circadian cycle
delW = .2     # depreciation rate per period for sleep capital when working
delP = .1     # depreciation rate per period for sleep capital when playng
gamma = .75   # curvature of consumption utility, CES
chiS = 5.     # utility weight on sleep timing
beta = .9     # subjective discount factor per period
# set up sine way for circadian cycle
yvect = np.linspace(0., 2*np.pi, num = S+1)
yvect = -np.sin(yvect + .5*b/np.pi)
plt.plot(yvect)
plt.show()
# save to mparams list
mparams = (S, eta, h, wbar, f, g, b, delW, delP, gamma, chiS, beta, yvect)

# set up other grids
nB = 11
nZ = 11
nC = 11
Bvect = np.linspace(0., 5., num = nB)
Zvect = np.linspace(0., 5., num = nZ)
Cvect = np.linspace(0., 5., num = nC)

# set up value function and jump functions
V = np.ones((nB,nZ,S))*(1/beta)
LamA = np.zeros((nB,nZ,S))
LamC = np.zeros((nB,nZ,S))
PhiB = np.zeros((nB,nZ,S))
PhiZ = np.zeros((nB,nZ,S))
Vnew = np.zeros((nB,nZ,S))

# set iteration parameters
maxiters = 500.
ccrit = 1.E-6

dist = 1.
iters = 0
while (dist > ccrit) and (iters < maxiters):
    iters = iters + 1
    for iy in range(0, S):
        for iB in range(0, nB):
            for iZ in range(0, nZ):
                maxval = -1.E+98
                for iA in range(0, 3):
                    for iC in range(0, nC):
                        # find utility, next period's savings and sleep stock
                        U, Bp, Zp = moddefs(Bvect[iB], Zvect[iZ], yvect[iy], \
                            iA, Cvect[iC], *mparams)
                        # find closest index for Bp and Zp
                        iBp, temp = nearest(Bvect, Bp)
                        iZp, temp = nearest(Zvect, Zp)
                        # next period's value for the circadian cycle
                        if iy == S-1:
                            iyp = 0
                        else:
                            iyp = iy + 1
                        # find new value
                        val = U + beta*V[iBp, iZp, iyp]
                        # if this is greater than previous max, update
                        if val > maxval:
                            maxval = val
                            Vnew[iB, iZ, iy] = val
                            LamA[iB, iZ, iy] = iA
                            LamC[iB, iZ, iy] = Cvect[iC]
                            PhiB[iB, iZ, iy] = Bvect[iBp]
                            PhiZ[iB, iZ, iy] = Zvect[iZp]
    dist = np.max(np.abs(V - Vnew))
    # replce old value function
    V = 1.*Vnew
    Vnew = 0.*Vnew
    print iters, dist             


# Simulate
T = 2*S    # number of periods to simulate

Bhist = np.zeros(T)
Zhist = np.zeros(T)
yhist = np.zeros(T)
Ahist = np.zeros(T)
Chist = np.zeros(T)

Bhist[0] = 1.0
Zhist[0] = 1.0
iy = 0

for t in range(0, T-1):
    # find closest index for B and Z
    iB, temp = nearest(Bvect, Bhist[t])
    iZ, temp = nearest(Zvect, Zhist[t])
    # find values for A and C
    Ahist[t] = LamA[iB, iZ, iy]
    Chist[t] = LamC[iB, iZ, iy]
    # update values of B, Z and y
    Bhist[t+1] = PhiB[iB, iZ, iy]
    Zhist[t+1] = PhiZ[iB, iZ, iy]
    if iy == S-1:
        iy = 0
    else:
        iy = iy + 1    
        
# plot data
t = range(0, T-1)
plt.plot(t, Zhist[0:T-1], label='Z')
plt.plot(t, Bhist[0:T-1], label='B')
plt.xlabel('time')
plt.show()

plt.plot(t, Ahist[0:T-1], label='A')
plt.xlabel('time')
plt.show()