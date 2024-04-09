import numpy as np

f= lambda x: x**2
a=0
b=5
n=5


def trap(f,n,a,b):
    x=np.linspace(a,b,n)
    y=f(x)
    total =y[0]+y[n-1]
    for j in range (n-2):
        total=total+2*y[j+1]
    total=total/(n-1)
    return total
    


trapz=trap(f,n,a,b)
        