# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 12:12:25 2017

@author: Kerk Phillips
"""

import numpy as np
import scipy.optimize as opt

def moddefs(B, Z, Bp, N, a, s, *mparams):
    '''
    B is the stock of savings today
    Z is the stock of sleep capital today
    Bp is the stock of savings tomorrow
    N is labor today
    a is wakeup time today
    s is bedtime today
    '''
    # labor effectiveness
    if eta == 1:
        e = np.log(Z)
    else:
        e  = h*(Z**(1-eta)-1)/(1-eta)
    # vertical shifter for circadian cycle
    d = np.arctan(f*Z)/np.pi + g
    # unconstrained wake and sleep times
    astar = np.pi + np.arcsin(d) - b
    sstar = 2*np.pi - np.arcsin(d) - b
    # consumption
    C = w*e*N + (1+r)*B - Bp
    # waking time
    T = (s - a)
    # leisure time
    L = T - N
    # marginal utility of consumption
    uC = C**(-gamma)
    # marginal utility of leisure
    uL = N**theta
    # marginal utilities of waking and sleeping
    ua = chis*(np.sin(a + b) - d)
    us = chis*(d - np.sin(s + b))
    # marginal utility of sleep stock
    uZ = chis*(astar - a + sstar - s)*f/(np.pi*(1+Z**2))
    # marginal effect of sleep stock on labor effectiveness
    phiZ = Z**(-eta)
    # utility from consumption
    if gamma == 1:
        UC = np.log(C)
    else:
        UC = (C**(1-gamma)-1)/(1-gamma)
    # utility from labor
    UL = - chiN*N**(1+theta)/(1+theta)
    # utiltiy from timing of sleep
    UT = chis*(d*(astar-a) + np.cos(astar+b) - np.cos(a-b) + \
               d*(sstar-s) + np.cos(s+b) - np.cos(sstar+b))
    # total utility
    U = UC + UL + UT
    return C, L, T, e, d, astar, sstar, uC, uL, ua, us, uZ, phiZ, U
    

def dynamic(invec, *mparams):
    # unpack invec
    Bpp = invec[0]  # savings in t+2
    Zpp = invec[1]  # sleep stock in t+2
    Bp = invec[2]   # savings is t+1
    Zp = invec[3]   # sleep stock in t+1
    B = invec[4]    # savings in t
    Z = invec[5]    # sleep stock in t
    Np = invec[6]   # labor in t+1
    ap = invec[7]   # wakeup time in t+1
    sp = invec[8]   # bedtime in t+1
    N = invec[9]    # labor in t
    a = invec[10]   # wakeup time in t
    s = invec[11]   # bedtime in t
    # get definitions for current period
    C, L, T, e, d, astar, sstar, uC, uL, ua, us, uZ, phiZ, U = \
        moddefs(B, Z, Bp, N, a, s, *mparams)
    # get definitions for next period
    Cp, Lp, Tp, ep, dp, astarp, sstarp, uCp, uLp, uap, usp, uZp, phiZp, Up = \
        moddefs(Bp, Zp, Bpp, Np, ap, sp, *mparams)  
    # labor-leisure Euler
    e1 = uL - uC*w*e
    # wake and sleep time Eulers
    e2 = ua + beta*(uCp*w*Np*phiZp + uZp)
    e3 = us - beta*(uCp*w*Np*phiZp + uZp)
    # intertemporal Euler
    e4 = uC - beta*uCp*(1+r)
    # sleep capital accumulation
    e5 = Zp - (1-deltaS)*Z + (1-T)
    return np.array([e1, e2, e3, e4, e5])


def steady(bar, *SSparams):
    # unpack bar
    [abar, sbar, Zbar, Cbar, astarbar, sstarbar, chiN, chis] = bar
    abar = abar % 24 # remainder after dividing by 24
    sbar = sbar % 24
    if sbar < abar:
        sbar = sbar + 24
    # definbe effectivness of labor
    ebar = h*(Zbar**(1-eta)-1)/(1-eta)
    # Euler equations and other definitions
    e1 = chiN*Nbar**theta - Cbar**(-gamma)*w*ebar
    temp = beta*(Cbar**(-gamma)*w*Nbar*Zbar**eta + chis*(astarbar - abar + \
        sstarbar - sbar)*f/(np.pi*(1+Zbar**2)))
    e2 = chis*(np.sin((abar + b)*np.pi/12) - dbar) + temp
    e3 = chis*(np.sin((abar + b)*np.pi/12) - dbar) - temp
    e4 = Zbar *deltaS - (24 - sbar + abar)
    e5 = Cbar - w*ebar*Nbar - r*Bbar
    e6 = dbar - np.arctan(f*Zbar)/np.pi + g
    e7 = astarbar - 12 + 12*np.sin(dbar)/np.pi + b
    e8 = sstarbar - 24 - 12*np.sin(dbar)/np.pi + b 
    return np.array([e1, e2, e3, e4, e5, e6, e7, e8])

    
# main program

# define parameters
eta = 2.      # curvature of phi(Z) functuion
h = 10.       # scaling factor for phi(Z) function
w = 1.        # constant wage
r = 0.        # constant interest rate per day
f = .5        # scaling factor for d(Z) function
g = .1        # additive factor for d(Z) function
b = 4.        # phase shifter for cyrcadian cycle
deltaS = .1   # depreciation rate per day for sleep capital
gamma = 2.5   # curvature of consumption utility, CES
chiN = .00014  # utility weight on leisure
chis = 60000. # utility weight on sleep timing
theta = .1    # curvature of labor disutility, CFE
beta = 1.     # subjective discount factor per day
Bbar = 0.
mparams = (eta, h, w, r, f, g, b, deltaS, gamma, chiN, chis, theta, beta)

# find steady state
# initial guesses
Nbar = 8
abar = 6
sbar = 22
Zbar = 88
Cbar = 80.
dbar = 0.42
astarbar = 7
sstarbar = 23.9

# setup things that are viewed as fixed paramters when finding the SS
SSparams = (eta, h, w, r, b, deltaS, gamma, chiN, chis, theta, beta, \
            Bbar, Nbar, dbar)
guess = np.array([abar, sbar, Zbar, Cbar, astarbar, sstarbar, chiN, chis])
# setup lambda function for use with fsolve
func = lambda bar: steady(bar, *SSparams)
# find the roots
bar = opt.fsolve(func, x0=guess)
print 'bar:  ', bar
check = steady(bar)
print 'check:', check

# unpack bar
[abar, sbar, Zbar, Cbar, astarbar, sstarbar, chiN, chis] = bar
abar = abar % 24  # divide by 24 return remainder
sbar = sbar % 24
if sbar < abar:
    sbar = sbar + 24
Tbar = sbar - abar
Lbar = Tbar - Nbar
Sbar = 24 - Tbar
print 'Nbar:', Nbar
print 'Lbar:', Lbar
print 'Sbar:', Sbar
print 'abar:', abar
print 'sbar:', sbar
print 'Zbar:', Zbar
print 'Cbar:', Cbar
print 'dbar:', dbar
print 'chiN:', chiN
print 'chis:', chis