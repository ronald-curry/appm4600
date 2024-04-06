import numpy as np
import matplotlib.pyplot as plt


def get_nodes(n,h):
    x=np.zeros(n+1)
    for i in range(n+1):
        temp1=(2*i+1)*np.pi
        temp2=2*(n+1)
        x[i]=np.cos(temp1/temp2)
    return x
def get_nodes_e(n,h):
    x=np.zeros(n+1)
    for i in range(n+1):
        x[i]=-1+(i-1)*h+h
    return x
   

def phi(z,x):
    p=1
    for i in range(n+1):
        p=p*(z-x[i])
    return p
def phi_e(z,x):
    p=1
    for i in range(n+1):
        p=p*(z-x[i])
    return p


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
            top=top+(w[j]/(z-x[j]))*f(x[j])
            bottom=bottom+(w[j]/(z-x[j]))
        if z==x[j]:
            err=j
    pz=top/bottom
    if err >-4:
        pz=f(x[err])
    return pz


f= lambda x: 1/(1+(16*x)**2)
n=80
h=2/n
x=get_nodes(n, h)
x_e=get_nodes_e(n, h)
#plt.plot(x,f(x),'o')
w=weights(n, x)
w_e=weights(n,x_e)
points=np.linspace(-1, 1,1001)
data=f(points)
poly=np.zeros(1001)
poly_e=np.zeros(1001)
ph=np.zeros(1001)
ph_e=np.zeros(1001)
for i in range(1001):
    poly[i]=p(x,f,w,n,points[i])
    poly_e[i]=p(x_e,f,w_e,n,points[i])
    ph[i]=phi(points[i],x)
    ph_e[i]=phi_e(points[i],x_e)
#plt.plot(points,data)
#plt.plot(points,poly)
#plt.plot(points,np.log10(abs(data-poly)))
#plt.plot(points,np.log10(abs(data-poly_e)))
#plt.show()

plt.plot(points,np.log10(ph))
plt.plot(points,np.log10(ph_e))
