import numpy as np


def eval_legendre(N,x):
    phi=np.zeros(N+1)
    phi[0]=1;
    phi[1]=x;
    for i in range(2,N+1):
        phi[i]=1/(i+1)*((2*i+1)*x*phi[i-1]-i*phi[i-1])
    p=phi;
    return p


legend=eval_legendre(6,2);