clc; close all; clear;


x0=1;
y0=1;
nMax=100;
tol=10^(-15);

guess1=fixedpoint(x0,y0,nMax,tol);

guess2=newton(x0,y0,nMax,tol);

guess3=Lazy(x0,y0,nMax,tol);
J(x0,y0)
function [z]=f(x,y)
    z=3*x.^2-y.^2;
end

function [z]=g(x,y)
    z=3*x.*y.^2-x.^3-1;
end

function [j]=J(x,y)
    j=[6*x,-2*y;3*y^2-3*x^2,6*x*y];
end

function guess=fixedpoint(x0,y0,nMax,tol)
    x=NaN(2,nMax+1);
    error=NaN(3,nMax+1);
    x(:,1)=[x0;y0];
    for n=1:nMax
        x(:,n+1)=x(:,n)-[1/6,1/18;0,1/6]*[f(x(1,n),x(2,n));g(x(1,n),x(2,n))];
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
    j = J(x(1,1),x(2,1));
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

