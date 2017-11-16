# -*- coding: utf-8 -*-
"""
This program is simulates a Two-Process model of sleep over a specified length
of time and plots the results.

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

def simulate(Sstart, Hstart, T, hrslate, *params):
    # given starting values, number of observations and parameters, calculate
    #  values of homeostatic cycle, sleep state, bounds, and circadian cycle
    # if hrslate is not 0, then stay up late on 2nd day
    if hrslate != 0:
        late = True
        latemin = ppday
        latemax = latemin + hrslate*ppday/24
    else:
        late = False
    
    yhist = np.zeros(T)
    for t in range(0, T):
        yhist[t] = circadian(t, *params)
        
    # calculate the upper and lower bounds
    Hphist = Hp + a*yhist
    Hmhist = Hm + a*yhist  
    
    # simulate
    Shist = np.zeros(T)
    Hhist = np.zeros(T)
    Shist[0] = Sstart
    Hhist[0] = Hstart
    for t in range(0, T-1):
        # determine state
        if late and (t >= latemin) and (t <= latemax) :
            Shist[t+1] = 1
        elif (Shist[t]==1) and (Hhist[t] < Hphist[t]):
            Shist[t+1] = 1
        elif (Shist[t]==1) and (Hhist[t] >= Hphist[t]):
           Shist[t+1] = 0 
        elif (Shist[t]==0) and (Hhist[t] > Hmhist[t]):
            Shist[t+1] = 0 
        else:
            Shist[t+1] = 1 
        # determine homestatic value
        if Shist[t+1]==1:
            Hhist[t+1] = mu + (Hhist[t] - mu)*np.exp(-1/chiW)
        else:
            Hhist[t+1] = Hhist[t]*np.exp(-1/chiS)
    return Hhist, Shist, Hphist, Hmhist, yhist  

def findSS(Sstart, Hstart, *params):
    # given starting values, output the value of H 24 hours later
    T = ppday+1
    Hhist, Shist, Hphist, Hmhist, yhist = \
       simulate(Sstart, Hstart, T, 0, *params)
    return Hhist[T-1]

# Main Program
# declare parameters
ndays = 3     # number of days to simulate
ppday= 24*60  # number of periods in a day
simple = 0    # if true use simple sine wave
alpha = 12.   # phase shift in hours
Hp = .60      # upper limit mean .60
Hm = .15      # lower limit mean .15
a = .1        # weight on circadian cycle .1
mu = 1.       # assmptote for upper limit 1.
chiS = 4.2    # half-life for sleep in hours 4.2
chiW = 18.2   # half-life for waking in hours 18.2
xhrs = 0      # if true show hours on x-axis, otherwise days
hrslate = 4.   # number of hours up late on 2nd day

# adjust and calculate paramters
omega = 2*np.pi/ppday  # ppday values per cycle rather than 2 pi
T = ppday*ndays + 1    # number of periods in simulation
chiS = chiS*ppday/24   # convert half life from hours to periods
chiW = chiW*ppday/24   # convert half life from hours to periods
alpha = alpha*ppday/24 # convert phase shift from hours to periods
if simple:
    avec = [1., 0., 0., 0., 0.]
else:
    avec = [.97, .22, .007, .03, .001] # coefficents on sine harmonics

# save parameters to a list
params = (alpha, Hp, Hm, a, mu, chiS, chiW, omega, avec)

# find steady state daily pattern, H the same 1 day later
dist = 1.
ccrit = 1.e-20
niter = 0
maxiter = 50
# starting values
Sstart = 0
Hstart = (Hp + Hm)/2
# apply contraction mapping
while (dist > ccrit) and (niter < maxiter) :
    niter += 1
    x = findSS(Sstart, Hstart, *params)
    dist = np.abs(Hstart - x)
    # print niter, dist, Hstart, x
    Hstart = x
print('Hstart is ', Hstart)
print('distance is ', dist)
print('number of iterations was ', niter)
# simulate
Hhist, Shist, Hphist, Hmhist, yhist = \
    simulate(Sstart, Hstart, T, 0, *params)

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
plt.plot(time, Hphist, label='upper bound')
plt.plot(time, Hmhist, label='lower bound')
plt.plot(time, Hhist, label='homestatic')
plt.legend(loc = 9, ncol = 3, fontsize = 10)
plt.ylim(0., .85)
plt.savefig('TwoProcess.eps', format='eps', dpi=2000)
plt.xlabel(xtext)
plt.show()

if hrslate !=0:
    # find new time path
    Hhist2, Shist2, temp1, temp2, temp3 = \
    simulate(Sstart, Hstart, T, hrslate, *params)
    
    # plot
    plt.figure()
    plt.plot(time, Hphist)
    plt.plot(time, Hmhist)
    plt.plot(time, Hhist, label='normal')
    plt.plot(time, Hhist2, label='up late')
    plt.legend(loc = 9, ncol = 2, fontsize = 10)
    plt.ylim(0., .85)
    plt.savefig('TwoProcess_late.eps', format='eps', dpi=2000)
    plt.xlabel(xtext)
    plt.show()