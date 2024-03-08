clc; close all; clear;


x0=[1.3;0;0];
nMax=30;
tol=10^(-6);
tol1=5*10^(-2);

tic;
newtong=newton(x0,nMax,tol);
toc

tic;
steepestg=steepest(x0,nMax,tol);
toc; 

tic;
broydeng=broyden(x0,nMax,tol);
toc;

tic;
steep_newtong=steepest_newton(x0,nMax,tol1,tol);
toc;

function [o]=f(x)
    o=x(1,:)+cos(x(1,:).*x(2,:).*x(3,:))-1;
end

function [o]=g(x)
    o=(1-x(1,:)).^(1/4)+x(2,:)+0.05*x(3,:).^2-0.15*x(3,:)-1;
end

function [o]=h(x)
    o=-x(1,:).^2-0.1*x(2,:).^2+0.01*x(2,:)+x(3,:)-1;
end

function [j]=J(x)
    j=[1-sin(x(1,:).*x(2,:).*x(3,:)).*x(2,:).*x(3,:),-sin(x(1,:).*x(2,:).*x(3,:)).*x(1,:).*x(3,:),-sin(x(1,:).*x(2,:).*x(3,:)).*x(2,:).*x(1,:);
    (1/4)*(1-x(1,:)).^(-3/4),1,0.1*x(3,:)-0.15;
    -2*x(1,:),-0.2*x(2,:)+0.01,1];
end

function [guess]=fixedpoint(x0,y0,nMax,tol)
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

function [guess]=newton(x0,nMax,tol)
    fun=@(x) ([f(x);g(x);h(x)]);
    Jfun=@(x) J(x);
    verb=0;
    [guess,x] = newton_method_nd(fun,Jfun,x0,tol,nMax,verb);
    x=cat(2,x0,x);
    figure()
    error=NaN(4,length(x(1,:)));
    error(1,:)=abs(g(x));
    error(2,:)=abs(g(x));
    error(3,:)=abs(h(x));
    error(4,:)=sqrt(error(1,:).^2+error(2,:).^2+error(3,:).^2);
    plot(1:length(x(1,:)),log10(error(4,:)))
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

function [guess]=broyden(x0,nMax,tol)
    B0=J(x0);
    fun=@(x) ([f(x);g(x);h(x)]);
    Bmat='fwd';
    verb=0;
    [guess,x] = broyden_method_nd(fun,B0,x0,tol,nMax,Bmat,verb);
    x=cat(2,x0,x);
    figure()
    error=NaN(4,length(x(1,:)));
    error(1,:)=abs(g(x));
    error(2,:)=abs(g(x));
    error(3,:)=abs(h(x));
    error(4,:)=sqrt(error(1,:).^2+error(2,:).^2+error(3,:).^2);
    plot(1:length(x(1,:)),log10(error(4,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Broyden")
end

function [guess2]=steepest_newton(x0,nMax,tol1,tol2)
    fun1=@(x) ([f(x);g(x);h(x)]);
    fun2=@(x) norm([f(x);g(x);h(x)]);
    Jfun=@(x) J(x);
    verb=0;
    Gfun=@(x) ([f(x);g(x);h(x)]);
    type="wolfe";
    [guess1,x1] = steepest_descent(fun2,Gfun,x0,tol1,nMax/2,type,verb);
    [guess2,x2] = newton_method_nd(fun1,Jfun,guess1,tol2,nMax/2,verb);
    x=cat(2,x1,x2);
    figure()
    error=NaN(4,length(x(1,:)));
    error(1,:)=abs(g(x));
    error(2,:)=abs(g(x));
    error(3,:)=abs(h(x));
    error(4,:)=sqrt(error(1,:).^2+error(2,:).^2+error(3,:).^2);
    plot(1:length(x(1,:)),log10(error(4,:)))
    xlabel("n")
    ylabel("log_{10}(error)")    
    title("Combined Newton/Steepest Descent")
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

function [guess]=steepest(x0,nMax,tol)
    fun=@(x) norm([f(x);g(x);h(x)]);
    Gfun=@(x) ([f(x);g(x);h(x)]);
    type="wolfe";
    verb=0;
    [guess,x] = steepest_descent(fun,Gfun,x0,tol,nMax,type,verb);
    x=cat(2,x0,x);
    figure()
    error=NaN(4,length(x(1,:)));
    error(1,:)=abs(g(x));
    error(2,:)=abs(g(x));
    error(3,:)=abs(h(x));
    error(4,:)=sqrt(error(1,:).^2+error(2,:).^2+error(3,:).^2);
    plot(1:length(x(1,:)),log10(error(4,:)))
    xlabel("n")
    ylabel("log_{10}(error)")
    title("Steepest Descent")
end

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
end