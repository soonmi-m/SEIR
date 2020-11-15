#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 19:16:41 2020

@author: soonmi
"""

import numpy as np
import matplotlib.pyplot as plt

N = 8880000 #population of nj
I0 = 98 #number of people infected at the start
E0 = 2100 #number of people exposed at the start
R0 = 2 #number of people recovered at the start
S0 = N - I0 - E0 - R0 #number of people suspectible at the start

Bet = .1 #beta - transmission rate
Gam = .04 #gamma - recovery rate
Sig = 0.20 #sigma - incubation rate
Mu = 0.000025 #mu - birth rate (increase in population)
Nu = 0.000025 #nu - death rate 

t = 40000 #total time 60 days

t0 = 1 #day 1

h = 1
n = (int)((t - t0)/h)  #number of iterations


def dSdt(t, S, E, I, R):
    return(Mu*N - Nu*S - Bet*S*I/N)

def dEdt(t, S, E, I, R):
    return(Bet*S*I/N - Nu*E -Sig*E)

def dIdt(t, S, E, I, R):
    return(Sig*E - Gam*I - Nu*I)

def dRdt(t, S, E, I, R):
    return(Gam*I - Nu*R)
    
S = S0
E = E0
I = I0
R = R0

t_list = np.zeros(n+1)
S_list = np.zeros(n+1)
E_list = np.zeros(n+1)
I_list = np.zeros(n+1)
R_list = np.zeros(n+1)
    
S_list[0] = S0
E_list[0] = E0
I_list[0] = I0
R_list[0] = R0


for i in range (1, n+1):
    s1 = h * dSdt(t0, S, E, I, R)
    s2 = h * dSdt(t0 + 0.5 * h, S + 0.5 * s1, E, I, R)
    s3 = h * dSdt(t0 + 0.5 * h, S + 0.5 * s2, E, I, R)
    s4 = h * dSdt(t0 + h, S + s3 , E, I, R)    
    
    
    e1 = h * dEdt(t0, S, E, I, R)
    e2 = h * dEdt(t0 + 0.5 * h, S, E + 0.5 * e1, I, R)
    e3 = h * dEdt(t0 + 0.5 * h, S, E + 0.5 * e2, I, R)
    e4 = h * dEdt(t0 + h, S, E + e3, I, R)

    
    i1 = h * dIdt(t0, S, E, I, R)
    i2 = h * dIdt(t0 + 0.5 * h, S, E, I + 0.5 * i1, R)
    i3 = h * dIdt(t0 + 0.5 * h, S, E, I + 0.5 * i2, R)
    i4 = h * dIdt(t0 + h, S, E, I + i3 , R)


    r1 = h * dRdt(t0, S, E, I, R)
    r2 = h * dRdt(t0 + 0.5 * h, S, E, I, R + 0.5 * r1)    
    r3 = h * dRdt(t0 + 0.5 * h, S, E, I, R + 0.5 * r2)
    r4 = h * dRdt(t0 + h, S , E, I, R + r3)    
    
    S = S + (1.0/6.0)*(s1 + 2 * s2 + 2 * s3 + s4)
    E = E + (1.0/6.0)*(e1 + 2 * e2 + 2 * e3 + e4)
    I = I + (1.0/6.0)*(i1 + 2 * i2 + 2 * i3 + i4)
    R = R + (1.0/6.0)*(r1 + 2 * r2 + 2 * r3 + r4)
   
    t_list[i] = i
    S_list[i] = S
    E_list[i] = E
    I_list[i] = I
    R_list[i] = R
    
    t0 = t0 + h


plt.plot(t_list, S_list, color = 'red', label= "Susceptible")
plt.plot(t_list, E_list, color = 'blue', label= "Exposed")
plt.plot(t_list, I_list, color = 'green', label= "Infected")
plt.plot(t_list, R_list, color = 'purple', label= "Recovered")

plt.title("New Jersey SEIR Model")
plt.xlabel("Days after March 15,2020")
plt.ylabel("Number of people")

plt.legend()