function [r,rn] = broyden_method_nd(fun,B0,x0,tol,nmax,Bmat,verb)
%{
This Matlab function finds, if possible, a root for the multivariable, 
differentiable function fun. Broyden method runs until
|fun(xn)|<tol or nmax iterations are reached. verb is a boolean (0 or 1)
controlling whether info is printed or not on the command line. 

Inputs: 
fun - (function or function handle) multivariable function fun(x), x is
assumed to be an n x 1 array
B0 - (n x n double) Initial "Jacobian-like" matrix B0 or its inverse  
x0   - (n x 1 double) initial guess for root
tol - (double) target accuracy / tolerance for the algorithm
nmax - (int) max number of iterations
Bmat - (string) Specifies one of three 'modes' for the first guess matrix
B0: 
        Bmat='fwd': B0 is an initial guess for the Jacobian Jf(x0)
        Bmat='inv': B0 is an initial guess for the inverse Jf(x0)^{-1}
        Bmat='eye' (or anything else): makes B0 the identity matrix
verb - (bool) print and pause (1) or don't print or pause (0)

Outputs:
r   - (n x 1 double) final estimate for the root
rn  - (n x niter double array) vector of iterates xn (used mainly for convergence /
testing purposes). 

Instructor: Eduardo Corona
%}

% verbose default is 0 
if nargin<7
    verb=0; 
end

% Initialize iterates and function value
xn=x0; rn(:,1)=x0; Fn = fun(xn); npn=1; 

% Depending on the mode Bmat, we set up functions to 'apply' the inverse
% and its transpose
if strcmp(Bmat,'fwd')
    I0 = @(x) B0\x; I0T = @(x) B0.'\x;  
elseif strcmp(Bmat,'inv')
   I0 = @(x) B0*x;  I0T = @(x) B0.'*x;
else
   %I0 and its transpose is the identity, so no function is needed.  
   I0 = []; I0T = []; 
end

% Start the arrays Un and Vn so that the update is I0*x + Un*(Vn'*x)
Un = zeros(length(x0),0); Vn = zeros(length(x0),0);

n=0; 
if verb
fprintf('\n|--n--|----xn----|---|f(xn)|---|')
end
while npn>tol && n<=nmax
    if verb
    fprintf('\n|--%d--|%1.7f|%1.7f|',n,norm(xn),norm(Fn));  
    %pause; 
    end
    
    % Broyden step xn-xn-1 = -Bk\f(xn)
    dn = -Inapp(I0,Un,Vn,Fn);
    % Update xn
    xn = xn + dn; 
    npn = norm(dn); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Typical formula updating Bn with Broyden and In with Sherman-Morrison
    %dfn = fun(xn)-fun(xn-dn);
    %rsn = dfn-Bn*dn; 
    % Update to the forward 'Jacobian like' matrix Bn
    %Bn = Bn + rsn*dn' / (dn'*dn);
    % Update to the inverse matrix (using Sherman-Morrison)
    %In = In + ((dn-In*dfn)/(dn'*(In*dfn)))*(dn'*In);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update In using only the previous In-1 (equivalent to the previous
    % formula for In) 
    Fn = fun(xn);
    un = Inapp(I0,Un,Vn,Fn); 
    cn = dn'*(dn+un);
    % The end goal is to add the rank 1 u*v' update as the next columns of
    % Vn and Un, as is done in, say, the eigendecomposition
    Vn = [Vn Inapp(I0T,Vn,Un,dn)];
    Un = [Un -(1/cn)*un];   
    
    n=n+1; 
    rn(:,n)=xn; 
end
    
r=xn;
    
if npn>tol
   fprintf('\n Broyden method failed to converge, n=%d, res=%e\n',nmax,norm(Fn)); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = Inapp(I0,Un,Vn,x)
%{
Function that applies I0*x + Un*Vn.'*x depending on a few cases for the
inputs
%}

if isempty(I0)
   if isempty(Un)
       y = x; 
   else
       y = x + Un*(Vn.'*x); 
   end
else
    if isempty(Un)
       y = I0(x);
   else
       y = I0(x) + Un*(Vn.'*x);
   end
end

end