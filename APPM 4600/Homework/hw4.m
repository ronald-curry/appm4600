clc; close all; clear;


x0=1;
y0=-1;
nMax=400;
tol=10^(-10);
tic;
guess1=newton(x0,y0,nMax,tol);
toc
tic;
guess2=Lazy(x0,y0,nMax,tol);

toc;
tic;
guess3=broyden(x0,y0,nMax,tol);
toc;

function [o]=f(x,y)
    o=x.^2+y.^2-4;
end

function [o]=g(x,y)
    o=exp(x)+y-1;
end


function [j]=J(x,y)
    j=[2*x,2*y;exp(x),1];
end

function guess=fixedpoint(x0,y0,nMax,tol)
    x=NaN(2,nMax+1);
    error=NaN(3,nMax+1);
    x(:,1)=[x0;y0];
    for n=1:nMax
        x(:,n+1)=x(:,n)-[1,0;0,1]*[f(x(1,n),x(2,n));g(x(1,n),x(2,n))]; %can change identiy to other matrix
        if sqrt((f(x(1,n+1),x(2,n+1)))^2+(g(x(1,n+1),x(2,n+1)))^2)<tol %checks if norm(f(n),g(n))<tol
            guess=[x(1,n+1);x(2,n+1)]; %final guess [x;y]
            break
        end
        if n==nMax
            error('did not converge to tol in Nmax')
        end
    end
    figure()
    error(1,:)=abs(f(x(1,:),x(2,:)));
    error(2,:)=abs(g(x(1,:),x(2,:)));
    error(3,:)=sqrt(error(1,:).^2+error(2,:).^2);
    plot(1:nMax+1,log10(error(3,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Fixed Point")
end

function [guess]=newton(x0,y0,nMax,tol)
    x=NaN(2,nMax+1);
    error=NaN(3,nMax+1);
    x(:,1)=[x0;y0];
    for n=1:nMax
        j = J(x(1,n),x(2,n));
        F = [f(x(1,n),x(2,n));g(x(1,n),x(2,n))];
        x(:,n+1)=x(:,n)-j\F;
        if sqrt((f(x(1,n+1),x(2,n+1)))^2+(g(x(1,n+1),x(2,n+1)))^2)<tol %checks if norm(f(n),g(n))<tol
            guess=[x(1,n+1);x(2,n+1)]; %final guess [x;y]
            break
        end
        if n==nMax
            error('did not converge to tol in Nmax')
        end
    end
    figure()
    error(1,:)=abs(f(x(1,:),x(2,:)));
    error(2,:)=abs(g(x(1,:),x(2,:)));
    error(3,:)=sqrt(error(1,:).^2+error(2,:).^2);
    plot(1:nMax+1,log10(error(3,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Newton")  
end

function [guess]=Lazy(x0,y0,nMax,tol)
    x=NaN(2,nMax+1);
    error=NaN(3,nMax+1);
    x(:,1)=[x0;y0];
    j = J(x0,y0);
    for n=1:nMax
        F = [f(x(1,n),x(2,n));g(x(1,n),x(2,n))];
        x(:,n+1)=x(:,n)-j\F;
        if sqrt((f(x(1,n+1),x(2,n+1)))^2+(g(x(1,n+1),x(2,n+1)))^2)<tol %checks if norm(f(n),g(n))<tol
            guess=[x(1,n+1);x(2,n+1)]; %final guess [x;y]
            break
        end
        if n==nMax
            error('did not converge to tol in Nmax')
        end
    end
    figure()
    error(1,:)=abs(f(x(1,:),x(2,:)));
    error(2,:)=abs(g(x(1,:),x(2,:)));
    error(3,:)=sqrt(error(1,:).^2+error(2,:).^2);
    plot(1:nMax+1,log10(error(3,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Lazy Newton")  
end

function [guess]=broyden(x0,y0,nMax,tol)
    B0=J(x0,y0);
    fun=@(x) ([f(x(1),x(2));g(x(1),x(2))]);
    x=[x0;y0];
    Bmat='fwd';
    verb=0;
    [guess,x] = broyden_method_nd(fun,B0,x,tol,nMax,Bmat,verb);
    x=cat(2,[x0;y0],x);
    figure()
    error=NaN(3,length(x(1,:)));
    error(1,:)=abs(f(x(1,:),x(2,:)));
    error(2,:)=abs(g(x(1,:),x(2,:)));
    error(3,:)=sqrt(error(1,:).^2+error(2,:).^2);
    plot(1:length(x(1,:)),log10(error(3,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Broyden")
end

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