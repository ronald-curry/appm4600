function [r,rn] = newton_method_nd(fun,Jfun,x0,tol,nmax,verb)
%{
This Matlab function finds, if possible, a root for the multivariable, 
differentiable function fun. Newton method runs until
|fun(xn)|<tol or nmax iterations are reached. verb is a boolean (0 or 1)
controlling whether info is printed or not on the command line. 

Inputs: 
fun - (function or function handle) multivariable function fun(x), x is
assumed to be an n x 1 array
Jfun - (function or function handle) Jacobian matrix of size n x n  
x0   - (n x 1 double) initial guess for root
tol - (double) target accuracy / tolerance for the algorithm
nmax - (int) max number of iterations
verb - (bool) print and pause (1) or don't print or pause (0)

Outputs:
r   - (n x 1 double) final estimate for the root
rn  - (n x niter double array) vector of iterates xn (used mainly for convergence /
testing purposes). 

Instructor: Eduardo Corona
%}

if nargin<6
    verb=0; 
end

% Initialize iteration and function value
xn=x0; rn(:,1)=x0; Fn = fun(xn); npn=1; 

n=0; 
if verb
fprintf('\n|--n--|----xn----|---|f(xn)|---|')
end
while npn>tol && n<=nmax
    Jn=Jfun(xn);  
    
    if verb
    fprintf('\n|--%d--|%1.7f|%1.7f|',n,norm(xn),norm(Fn));  
    %pause; 
    end
    
    % Newton step x_{n+1} = x_n - Jf(x_n)^{-1} * F(x_n) 
    pn = -Jn\Fn;
    xn = xn + pn;
    npn = norm(pn); 
    
    n=n+1; 
    rn(:,n)=xn;
    Fn = fun(xn); 
end
    
r=xn;
    
if npn>tol
   fprintf('\n Newton method failed to converge, n=%d, res=%e\n',nmax,norm(Fn)); 
end

end