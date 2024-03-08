function [alpha,x1] = line_search(fun,Gfun,x0,p,type,params)
%{
This Matlab function performs a back-tracking line search to find a value 
of alpha and x1 = x0 + alpha*p such that sufficient descent conditions have
been achieved for function fun in the direction p.  

Inputs: 
fun - (function or function handle) function fun(x), x is
assumed to be an n x 1 array
Gfun - (function or function handle) Gradient vector of size n x 1
x0   - (n x 1 double) initial guess for minimum
p    - (n x 1 double) descent direction p 
type - (string) line search type. Options are 'wolfe', 'swolfe' (symmetric wolfe) and
'armijo' (default). 
params - (struct) parameters struct containing arguments for constants c1
and c2 (for descent and curvature conditions) and maxback (maximum number
of backtracking steps)

Outputs:
r   - (n x 1 double) final estimate for the minimum
rn  - (n x niter double array) vector of iterates xn (used mainly for convergence /
testing purposes). 

Instructor: Eduardo Corona
%}

alpha=2; n=0; nmax = params.maxback; 
cond = false; 
f0 = fun(x0); 
Gdotp = Gfun(x0)'*p;

% We run a while loop until the required descent conditions, encoded in the
% boolean variable cond are met, or n reaches nmax (maximum number of
% backtracking steps)
while n<=nmax & cond==false
    alpha = alpha/2; 
    x1 = x0+alpha*p;
if strcmp(type,'wolfe')
    % Wolfe (Armijo sufficient descent and simple curvature conditions)
    Armijo = fun(x1) <= f0 + params.c1*alpha*Gdotp; 
    Curvature = Gfun(x1)'*p >= params.c2*Gdotp;
    cond = Armijo & Curvature; 
elseif strcmp(type,'swolfe')
    % Symmetric Wolfe (Armijo sufficient descent and symmetric curvature conditions)
    Armijo = fun(x1) <= f0 + params.c1*alpha*Gdotp; 
    Curvature = abs(Gfun(x1)'*p) <= params.c2*abs(Gdotp);
    cond = Armijo & Curvature;
else
    % Armijo (only sufficient descent condition)
    Armijo = fun(x1) <= f0 + params.c1*alpha*Gdotp;
    cond = Armijo; 
end

n=n+1; 

end