import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import scipy.linalg as scila
import time

def driver():

     ''' create  matrix for testing different ways of solving a square 
     linear system'''

     '''' N = size of system'''
     N = 1000
 
     ''' Right hand side'''
     '''
     b = np.random.rand(N,1)
     A = np.random.rand(N,N)
     t1=time.time()
     x = scila.solve(A,b)
     t2=time.time()
     test = np.matmul(A,x)
     r = la.norm(test-b)
     print("Regular solve")
     print(r)
     print(t2-t1)
     t3=time.time()
     LU,P=scila.lu_factor(A)
     t4=time.time()
     x2 = scila.lu_solve((LU,P), b)
     t5=time.time()
     test2=np.matmul(A,x2)
     r2=la.norm(test2-b)
     print("LU")
     print(r2)
     print(t4-t3)
     print(t5-t4)
     '''
     ''' Create an ill-conditioned rectangular matrix '''
     N = 10
     M = 5
     A = create_rect(N,M)     
     b = np.random.rand(N,1)
     x=scila.solve(A.T@A,A.T@b)
     test = np.matmul(A,x)
     r = la.norm(test-b)
     print('Normal')
     print(r)
     Q,R=la.qr(A)
     x2=scila.solve(R,Q.T@b)
     test2=np.matmul(A,x2)
     r2=la.norm(test2-b)
     print('QR')
     print(r2)
     
     
def create_rect(N,M):
     ''' this subroutine creates an ill-conditioned rectangular matrix'''
     a = np.linspace(1,15,M)
     d = 10**(-a)
     
     D2 = np.zeros((N,M))
     for j in range(0,M):
        D2[j,j] = d[j]
     
     '''' create matrices needed to manufacture the low rank matrix'''
     A = np.random.rand(N,N)
     Q1, R = la.qr(A)
     test = np.matmul(Q1,R)
     A =    np.random.rand(M,M)
     Q2,R = la.qr(A)
     test = np.matmul(Q2,R)
     
     B = np.matmul(Q1,D2)
     B = np.matmul(B,Q2)
     return B     
          
  
if __name__ == '__main__':
      # run the drivers only if this is called from the command line
      driver()       
