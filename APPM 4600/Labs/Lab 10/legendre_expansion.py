import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad as quad

def driver():

#  function you want to approximate
    #f = lambda x: math.exp(x)
    f= lambda x: 1/(1+x**2)
# Interval of interest
    a = -1
    b = 1
# weight function
    w =  lambda x: 1
    w2 = lambda x: 1/np.sqrt(1-x**2)
# order of approximation
    n = 3

#  Number of points you want to sample in [a,b]
    N = 1000
    xeval = np.linspace(a,b,N+1)
    pval = np.zeros(N+1)

    for kk in range(N+1):
      pval[kk] = eval_legendre_expansion(f,a,b,w,n,xeval[kk])

    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])

    plt.figure();
    plt.plot(xeval,pval);
    plt.show()

    plt.figure();
    err = abs(pval-fex)
    plt.plot(xeval,np.log10(err)); 
    plt.show()
    
    pval = np.zeros(N+1)
    for kk in range(N+1):
      pval[kk] = eval_cheby_expansion(f,a,b,w2,n,xeval[kk])

    ''' create vector with exact values'''
    fex = np.zeros(N+1)
    for kk in range(N+1):
        fex[kk] = f(xeval[kk])

    plt.figure();
    plt.plot(xeval,pval);
    plt.show()

    plt.figure();
    err = abs(pval-fex)
    plt.plot(xeval,np.log10(err)); 
    plt.show()




def eval_legendre_expansion(f,a,b,w,n,x):

#   This subroutine evaluates the Legendre expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab
  p = eval_legendre(n, x)
  # initialize the sum to 0
  pval = 0.0
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_legendre(j, x)[j]
      func1 = lambda x: phi_j(x)*f(x)*w(x)
      func2 =lambda x: phi_j(x)**2*w(x)
      # use the quad function from scipy to evaluate coeffs
      funcj,err=quad(func1,a,b)
      norm,err=quad(func2,a,b)
      aj=funcj/norm
      # accumulate into pval
      pval = pval+aj*p[j]

  return pval

def eval_cheby_expansion(f,a,b,w,n,x):

#   This subroutine evaluates the Chebyshev expansion

#  Evaluate all the Legendre polynomials at x that are needed
# by calling your code from prelab
  p = eval_cheby(n, x)
  # initialize the sum to 0
  pval = 0.0
  for j in range(0,n+1):
      # make a function handle for evaluating phi_j(x)
      phi_j = lambda x: eval_cheby(j, x)[j]
      func1 = lambda x: phi_j(x)*f(x)*w(x)
      func2 =lambda x: phi_j(x)**2*w(x)
      # use the quad function from scipy to evaluate coeffs
      funcj,err=quad(func1,a,b)
      norm,err=quad(func2,a,b)
      aj=funcj/norm
      # accumulate into pval
      pval = pval+aj*p[j]

  return pval

def eval_legendre(N,x):
    phi=np.zeros(N+1)
    phi[0]=1;
    if N>0:
        phi[1]=x;
    if N>1:
        for i in range(2,N+1):
            phi[i]=1/(i)*((2*(i-1)+1)*x*phi[i-1]-(i-1)*phi[i-2])
    return phi

def eval_cheby(n,x):
    T=np.zeros(n+1)
    T[0]=1;
    if n>0:
        T[1]=x;
    if n>1:
        for i in range(2,n+1):
            T[i]=2*(i-1)*T[i-1]-T[i-2]
    return T

if __name__ == '__main__':
  # run the drivers only if this is called from the command line
  driver()
