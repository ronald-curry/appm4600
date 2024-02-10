# import libraries
import numpy as np
    
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
f1 = lambda x: x-(x**5-7)/(5*x**4)

Nmax = 11
tol = 1e-15

''' test f1 '''
x0 = 1
[xstar,ier,x] = fixedpt(f1,x0,tol,Nmax)
print('the approximate fixed point is:',xstar)
print('f1(xstar):',f1(xstar))
print('Error message reads:',ier)
    
