# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 12:14:34 2024

@author: ronal
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def fixedpt(f,x0,tol,Nmax):

    ''' x0 = initial guess''' 
    ''' Nmax = max number of iterations'''
    ''' tol = stopping tolerance''' 
    
    
    ''' Outputs:
    xstar: final guess of fixed point
    ier: error message
    x: array of all guesses
    '''
    x = np.zeros((Nmax+1,1))
    count = 0
    x[0]=x0
    while (count <Nmax):
        x1 = f(x0)
        count = count +1
        x[count]=x1
        if (abs(x1-x0) <tol):
            xstar = x1
            ier = 0
            return [xstar,ier,x,count]
        x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,x,count]

def bisection(f,a,b,tol,Nmax):
    '''
    Inputs:
      f,a,b       - function and endpoints of initial interval
      tol, Nmax   - bisection stops when interval length < tol
                  - or if Nmax iterations have occured
    Returns:
      astar - approximation of root
      ier   - error message
            - ier = 1 => cannot tell if there is a root in the interval
            - ier = 0 == success
            - ier = 2 => ran out of iterations
            - ier = 3 => other error ==== You can explain
    '''

    '''     first verify there is a root we can find in the interval '''
    fa = f(a); fb = f(b);
    count = 0
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier,count]

    ''' verify end point is not a root '''
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier,count]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier,count]

    
    while (count < Nmax):
      c = 0.5*(a+b)
      fc = f(c)

      if (fc ==0):
        astar = c
        ier = 0
        return [astar, ier,count]

      if (fa*fc<0):
         b = c
      elif (fb*fc<0):
        a = c
        fa = fc
      else:
        astar = c
        ier = 3
        return [astar, ier,count]

      if (abs(b-a)<tol):
        astar = a
        ier =0
        return [astar, ier,count]
      
      count = count +1

    astar = a
    ier = 2
    return [astar,ier,count] 

# use routines    
f = lambda x: 1
a = 1
b = 4

Nmax = 130
tol = 1e-16

[astar,ier,count] = bisection(f,a,b,tol,Nmax)
#print('the approximate root is',astar)
#print('the error message reads:',ier)
#print("number of iterations",count)

f2 = lambda x: 12/(x+1)
x0 = 2.8
[xstar,ier,x,count] = fixedpt(f2,x0,tol,Nmax)

print('the approximate fixed point is:',xstar)
print('f1(xstar):',f2(xstar))
print('Error message reads:',ier)
print("Number of iterations", count)
#print('Iterations',x)