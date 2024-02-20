# import libraries
import numpy as np

def in_basin (f,fp,fpp,x):
    y=f(x)*fpp(x)/(fp(x)**2)
    indicator=0
    if abs(y)<1:
        indicator=1
    return indicator
# define routines
def bisection(f,fp,fpp,a,b,tol,in_basin):
    
#    Inputs:
#     f,a,b       - function and endpoints of initial interval
#      tol  - bisection stops when interval length < tol

#    Returns:
#      astar - approximation of root
#      ier   - error message
#            - ier = 1 => Failed
#            - ier = 0 == success

#     first verify there is a root we can find in the interval 

    fa = f(a)
    fb = f(b);
    if (fa*fb>0):
       ier = 1
       astar = a
       return [astar, ier]

#   verify end points are not a root 
    if (fa == 0):
      astar = a
      ier =0
      return [astar, ier]

    if (fb ==0):
      astar = b
      ier = 0
      return [astar, ier]

    count = 0
    d = 0.5*(a+b)
    while (in_basin(f,fp,fpp,d)==0):
      fd = f(d)
      if (fd ==0):
        astar = d
        ier = 0
        return [astar, ier]
      if (fa*fd<0):
         b = d
      else: 
        a = d
        fa = fd
      d = 0.5*(a+b)
      count = count +1
#      print('abs(d-a) = ', abs(d-a))
      
    astar = d
    ier = 0
    return [astar, ier,count]


def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]
        
# use routine
f = lambda x: np.exp(x**2+7*x-30)-1
fp = lambda x: (2*x+7)*np.exp(x**2+7*x-30)
fpp = lambda x: (2)**np.exp(x**2+7*x-30)+(2*x+7)**2*np.exp(x**2+7*x-30)


Nmax = 100
tol=1e-10
a = 2
b = 4.5
[astar,ier,count] = bisection(f,fp,fpp,a,b,tol,in_basin)
print(astar)
p0 = astar
[p,pstar,info,it] = newton(f,fp,p0,tol, Nmax)
print(pstar)
print("Number of iterations",count+it)



