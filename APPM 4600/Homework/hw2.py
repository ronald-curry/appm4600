# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 19:18:53 2024

@author: ronal
"""

import numpy as np
import matplotlib.pyplot as plt
import math
'''
t=np.zeros(31);
#print(t[0])
for i in range(31):
    t[i]=i*np.pi/30
y=np.cos(t)
print("The sum, S, is:")
print(np.sum(y))
'''
'''
theta=np.linspace(0,2*np.pi)
R=1.2
delta_r=0.1
f=15
p=0
x=R*(1+delta_r*np.sin(f*theta+p))*np.cos(theta)
y=R*(1+delta_r*np.sin(f*theta+p))*np.sin(theta)
plt.plot(x,y)
plt.show()
for i in range (10):
    R=i
    delta_r=0.05
    f=2+1
    p=np.random.uniform(0,2)
    x=R*(1+delta_r*np.sin(f*theta+p))*np.cos(theta)
    y=R*(1+delta_r*np.sin(f*theta+p))*np.sin(theta)
    plt.plot(x,y)
plt.show()
'''

x=9.999999995000000 * 10**(-10)
y=math.e**x
y=y-1
print("y calc using e^")
print(y)
