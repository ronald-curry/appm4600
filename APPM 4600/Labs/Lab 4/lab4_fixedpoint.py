# import libraries
import numpy as np
    

def delta_2(f,slow,Nmax,tol):
    count = 0
    fast=np.zeros((Nmax,1))
    while (count < Nmax):
        fast[count]=slow[count]-(slow[count+1]-slow[count])**2/(slow[count+2]-2*slow[count+1]+slow[count])
        count=count+1
        if slow[count+2]==0:
            break
        if (f(fast[count])-fast[count]<tol):
            break
    return fast

# define routines
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
            return [xstar,ier,x]
        x0 = x1

    xstar = x1
    ier = 1
    return [xstar, ier,x]
    

# use routines 
f1 = lambda x: (10/(x+4))**0.5

Nmax = 20
tol = 1e-5

''' test f1 '''
x0 = 1.5
[xstar,ier,x] = fixedpt(f1,x0,tol,Nmax)
aitkens=delta_2(f1,x,Nmax,tol)

print('the approximate fixed point is:',xstar)
print('f1(xstar):',f1(xstar))
print('Error message reads:',ier)
print('Iterations',x)
print('Aitkens: ',aitkens)
    
