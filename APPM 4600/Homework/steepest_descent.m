function [r,rn] = steepest_descent(fun,Gfun,x0,tol,nmax,type,verb)
%{
This Matlab function finds, if possible, a local minimum for the  
differentiable function fun with gradient Gfun using the basic 
Gradient descent step starting from initial guess x0. 
Step size is restricted with a simple back-tracking line search strategy. 

It is generally assumed that the function is locally convex at x0 
(the Hessian Hfun(x0) is symmetric, positive definite) 
 
The method runs until |fun(xn)|<tol or nmax iterations are reached. 
type is a string indicating line search type and verb is a boolean (0 or 1)
controlling whether info is printed or not on the command line. 

Inputs: 
fun - (function or function handle) function fun(x), x is
assumed to be an n x 1 array
Gfun - (function or function handle) Gradient vector of size n x 1
x0   - (n x 1 double) initial guess for minimum
tol - (double) target accuracy / tolerance for the algorithm
nmax - (int) max number of iterations
type - (string) line search type. Options are 'wolfe', 'swolfe' and
'armijo' (default). 
verb - (bool) print and pause (1) or don't print or pause (0)

Outputs:
r   - (n x 1 double) final estimate for the minimum
rn  - (n x niter double array) vector of iterates xn (used mainly for convergence /
testing purposes). 

Instructor: Eduardo Corona
%}

params.c1 = 10^-3; params.c2 = 0.9; params.maxback=10; xn=x0; 
n=0; rn(:,1)=x0; 
alpha=1; 
fn = fun(xn); 
pn = -Gfun(xn); 

if verb
fprintf('\n|--n--|-alpha-|----|xn|----|---|f(xn)|---|---|Gf(xn)|---|')
end

while n<=nmax && norm(pn)>tol
    if verb
    fprintf('\n|--%d--|%1.5f|%1.7f|%1.7f|%1.7f|',n,alpha,norm(xn),abs(fn),norm(pn));  
    %pause(0.01); 
    end
    
    [alpha,xn] = line_search(fun,Gfun,xn,pn,type,params);
    fn = fun(xn); 
    pn = -Gfun(xn);
    n=n+1; 
    rn(:,n+1)=xn; 
end

r=xn; 

end