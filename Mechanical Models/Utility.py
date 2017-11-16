# -*- coding: utf-8 -*-
"""
This program is simulates a utility model of sleep over a specified length of 
time and plots the results.

Code written by Kerk l. Phillips
Mar. 11, 2017
"""

import numpy as np
import matplotlib.pyplot as plt

def cterm(t, i, *params):
    # caclulate the ith harmonic
    x = avec[i] * np.sin((i+1)*omega*(t-alpha))
    return x

def circadian(t, *params):
    # calculate and add up all 5 harmonics
    y = 0.
    for i in range(0, 5):
        y += cterm(t, i, *params)
    return y    

def simulate(Sstart, Zstart, T, hrslate, stime, *params):
    # given starting values, number of observations and parameters, calculate
    #  values of homeostatic cycle, sleep state, bounds, and circadian cycle
    # if hrslate is not 0, then stay up late on 2nd day
    if hrslate != 0:
        late = True
        latemin = float(ppday + stime)
        latemax = latemin + hrslate*ppday/24
    else:
        late = False
    
    yhist = np.zeros(T)
    for t in range(0, T):
        yhist[t] = circadian(t, *params)
        
    # simulate
    Zhist = np.zeros(T)
    Shist = np.zeros(T)
    dhist = np.zeros(T)
    Shist[0] = Sstart
    Zhist[0] = Zstart
    for t in range(0, T):
        dhist[t] = f*Zhist[t] + g
        # determine state
        if late and (t >= latemin) and (t <= latemax) :
            Shist[t] = 1
        elif yhist[t] + dhist[t] > 0.:
            Shist[t] = 1 
        else:
            Shist[t] = 0 
        # determine sleep stock value
        if t < T-1:
            if Shist[t] == 1:
                Zhist[t+1] = Zhist[t]*(1-delta)
            else:
                Zhist[t+1] = Zhist[t]*(1-delta) + 1./ppday
    
    return dhist, Zhist, Shist, yhist  

def findSS(Sstart, Zstart, *params):
    # given starting values, output the value of H 24 hours later
    T = ppday+1
    dhist, Zhist, Shist, yhist = \
       simulate(Sstart, Zstart, T, 0, *params)
    return Zhist[T-1]

def findnth(n, val, array, comp):
    # find the index where array execeeds val for the nth time
    if comp == 'gt':
        cond = array > val
    else:
        cond = array < val
    counts = np.cumsum(cond)
    idx = np.searchsorted(counts, n)
    return idx

def findsw(array):
    T = array.size
    diff = np.zeros(T-1)
    for t in range(0, T-1):
        diff[t] = array[t+1] - array[t]
    stime = findnth(1, 0, diff, 'lt')
    wtime = findnth(1, 0, diff, 'gt')
    return stime, wtime


# Main Program
# declare parameters
ndays = 3     # number of days to simulate
ppday= 24*60  # number of periods in a day
simple = 0    # if true use simple sine wave
alpha = 12.   # phase shift in hours
xhrs = 0      # if true show hours on x-axis, otherwise days
hrslate = 4.  # number of hours up late on 2nd day
delta = .05   # sleep stock depreciation rate per hour
f = 4.
g = -.3

# normal is delta - .05, f= 4., g = -.3
# fuzzy is delta = .1, f = 10., g = -1.45

# adjust and calculate paramters
omega = 2*np.pi/ppday  # ppday values per cycle rather than 2 pi
T = ppday*ndays + 1    # number of periods in simulation
alpha = alpha*ppday/24 # convert phase shift from hours to periods
delta = 1. - (1. - delta)**(24.0/ppday)
if simple:
    avec = [1., 0., 0., 0., 0.]
else:
    avec = [.97, .22, .007, .03, .001] # coefficents on sine harmonics

# save parameters to a list
params = (alpha, omega, f, g, delta, avec)

# find steady state daily pattern, H the same 1 day later
# parameter for the contraction mapping
dist = 1.
ccrit = 1.e-20
niter = 0
maxiter = 50
# starting values
Sstart = 0.
Zstart = 1.
# apply contraction mapping
while (dist > ccrit) and (niter < maxiter) :
    niter += 1
    Znew = findSS(Sstart, Zstart, *params)
    dist = np.abs(Zstart - Znew)
    print(niter, dist, Zstart, Znew)
    Zstart = Znew
print('Zstart is ', Zstart)
print('distance is ', dist)
print('number of iterations was ', niter)

# simulate
dhist, Zhist, Shist, yhist = \
    simulate(Sstart, Zstart, T, 0, 0, *params)
mdhist = -1.*dhist

# find normal sleep time
stime, wtime = findsw(Shist)

# plot cycle
# get time value
time = np.linspace(0, T, T+1)
# convert to hours or days
if xhrs:
    time = 24*time[0:T]/ppday
    xtext = 'time in hours'
else:
    time = time[0:T]/ppday
    xtext = 'time in days'
    
plt.figure()
plt.plot(time, mdhist, label='-d')
plt.plot(time, yhist, label='circadian')
plt.legend(loc = 9, ncol = 3, fontsize = 10)
# plt.ylim(0., .85)
plt.savefig('Utility.eps', format='eps', dpi=2000)
plt.xlabel(xtext)
plt.show()

plt.figure()
plt.plot(time, Zhist)
plt.savefig('Utility_Z.eps', format='eps', dpi=2000)
plt.xlabel(xtext)
plt.show()


if hrslate !=0:
    # find new time path
    dhist2, Zhist2, Shist2, temp2 = \
    simulate(Sstart, Zstart, T, hrslate, stime, *params)
    mdhist2 = -1.*dhist2
    
    # plot
    plt.figure()
    plt.plot(time, yhist)
    plt.plot(time, mdhist, label='normal')
    plt.plot(time, mdhist2, label='up late')
    plt.legend(loc = 9, ncol = 2, fontsize = 10)
    # plt.ylim(0., .85)
    plt.savefig('Utility_late.eps', format='eps', dpi=2000)
    plt.xlabel(xtext)
    plt.show()