# This script tests the convergence of adaptive quad
# and compares to a non adaptive routine
import numpy as np
# get adaptive_quad routine and numpy from adaptive_quad.py
from adaptive_quad import *
# get plot routines
import matplotlib.pyplot as plt
import numpy as np
# specify the quadrature method
# (eval_gauss_quad, eval_composite_trap, eval_composite_simpsons)

# interval of integration [a,b]
a = 0.1; b = 2.
# function to integrate and true values
# TRYME: uncomment and comment to try different funcs
#        make sure to adjust I_true values if using different interval!
#f = lambda x: np.log(x+(x==0))**2; I_true = 2; labl = '$\log^2(x)$'
#f = lambda x: 1./(np.power(x,(1./5.))); I_true = 5./4.; labl = '$\\frac{1}{x^{1/5}}$'
#f = lambda x: np.exp(np.cos(x)); I_true = 2.3415748417130531; labl = '$\exp(\cos(x))$'
#f = lambda x: x**20; I_true = 1./21.; labl = '$x^{20}$'
# below is for a=0.1, b = 2
#a=0.1; b=2;
f = lambda x: np.sin(1./x); I_true = 1.1455808341; labl = '$\sin(1/x)$'
# absolute tolerance for adaptive quad
tol = 1e-3
# machine eps in numpy
eps = np.finfo(float).eps




M=5
method = eval_gauss_quad
I_gauss,X_gauss,nsplit_gauss = adaptive_quad(a,b,f,tol,M,method)
method = eval_composite_trap
I_trap,X_trap,nsplit_trap = adaptive_quad(a,b,f,tol,M,method)
method = eval_composite_simpsons
I_simp,X_simp,nsplit_simp = adaptive_quad(a,b,f,tol,M,method)


