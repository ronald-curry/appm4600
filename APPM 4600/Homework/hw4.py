import numpy as np
import matplotlib.pyplot as plt
import math




def driver():

    x0=1
    y0=1
    nMax=20
    tol=10**(-5)
    
    x=np.array([[x0],[y0]])
    for n in range(nMax):
        x=x-np.array([[1/6,1/18],[0,1/6]])*[[f(x[0],x[1])],[g(x[0],x[1])]]
        if np.sqrt(x[0]**2+x[1]**2)<tol:
            break
        
    print(x)

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)
  
    


''' create divided difference matrix'''
def dividedDiffTable(x, y, n):
 
    for i in range(1, n):
        for j in range(n - i):
            y[j][i] = ((y[j][i - 1] - y[j + 1][i - 1]) /
                                     (x[j] - x[i + j]));
    return y;
    
def evalDDpoly(xval, xint,y,N):
    ''' evaluate the polynomial terms'''
    ptmp = np.zeros(N+1)
    
    ptmp[0] = 1.
    for j in range(N):
      ptmp[j+1] = ptmp[j]*(xval-xint[j])
     
    '''evaluate the divided difference polynomial'''
    yeval = 0.
    for j in range(N+1):
       yeval = yeval + y[0][j]*ptmp[j]  

    return yeval

       
def vandermonde(x,y,N):
    v=np.zeros((N+1,N+1))
    for i in range(0,N+1):
        for j in range(0,N+1):
            v[i,j]=x[i]**j
    a=np.linalg.solve(v, y)
    return a
  

def f(x,y):
    return 3*x**2-y**2

def g(x,y):
    return 3*x*y**2-x**3-1

driver()