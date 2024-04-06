import numpy as np
import matplotlib.pyplot as plt


def get_nodes(n,h):
    x=np.zeros(n+1)
    for i in range(n+1):
        x[i]=-1+(i-1)*h+2/n
    return x
   

def weights(n,x):
    w=np.zeros(n+1)
    for j in range(n+1):
        p=1
        for i in range(n+1):
            if j!=i:
                p=p*(x[j]-x[i])
        w[j]=1/p
    return w
def p(x,f,w,n,z):
    top=0
    bottom=0
    err=-5
    for j in range(n+1):
        if z!=x[j]:
            top=top+w[j]/(z-x[j])*f(x[j])
            bottom=bottom+w[j]/(z-x[j])
        if z==x[j]:
            err=j
    pz=top/bottom
    if err >-4:
        pz=f(x[err])
    return pz


f= lambda x: 1/(1+(16*x)**2)
n=5
h=2/n
x=get_nodes(n, h)
plt.plot(x,f(x),'o')
w=weights(n, x)
points=np.linspace(-1, 1,1001)
data=f(points)
poly=np.zeros(1001)
for i in range(1001):
    poly[i]=p(x,f,w,n,data[i])
plt.plot(points,data)
plt.plot(points,poly,'o')
plt.show()