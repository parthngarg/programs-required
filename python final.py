# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 02:20:22 2016

@author: rajivgarg
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 23:37:21 2016

@author: rajivgarg
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
from numpy import pi
a=10
c=1
rho=1
fi=1
ti=29
cp=1
i=10000
t = np.linspace(0,400,i)
q=rho*cp*scipy.sin(2*pi*t)


"function to plot" 
def stirredtank(ht,t):

    dhbydt=(fi-c*math.sqrt(ht[0]))/a
    
    q=rho*cp*fi*0.5*scipy.sin(0.1*t)    

    dt2bydt=(fi*(ti-ht[1]) + q/rho*cp)/(a*ht[0])

    ans=[dhbydt,dt2bydt]

    return ans
    
h=[2,18]
q=1
y=odeint(stirredtank,h,t)
print y[:,1]
    
a2=1.27*10**-4
rho2=1.355*10**4
v2=8.84*10**-8
cp2=139.4
h2=750
tau=rho2*v2*cp2/(h2*a2)
tsp=30

print "Instrumentation-control valves"
print "Linear diaphragm characteristic valves"
print "piston type pneumaatic actuators"
print "type E type thermocoupls- can be used in oxidizing and inert atmosphere"

v=10
j=0

tttt=10*scipy.sin(2*pi*t)
temp2=[]
yy=[]

while(j<i):
    yy.append(y[j,1])
    temp2.append(yy[j]*(1-np.exp(-t[j]/tau)))
    j=j+1

from sympy.integrals.transforms import inverse_laplace_transform
from sympy.integrals import laplace_transform
#from sympy.abc import fi, s,cp, tsp, rho, v/fi,tttt

#laplace_transform(t**a, t, s)

fig = plt.figure()
plt.plot(t,y[:,0],'b--')
fig.suptitle('height vs time', fontsize=20)
plt.xlabel('time', fontsize=18)
plt.ylabel('height', fontsize=16)

fig = plt.figure()
plt.plot(t,y[:,1],'r--')
fig.suptitle('temp vs time', fontsize=20)
plt.xlabel('time', fontsize=18)
plt.ylabel('temperature', fontsize=16)

fig = plt.figure()
plt.plot(t,temp2,'r--')
fig.suptitle('temp(as measured by thermometre) vs time', fontsize=20)
plt.xlabel('time', fontsize=18)
plt.ylabel('temperature', fontsize=16)
i=100

#if sinosoidal input is given to q