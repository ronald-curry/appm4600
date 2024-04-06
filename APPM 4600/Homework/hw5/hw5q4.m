clear; clc; close all;
x=[0 1 2 3];

M=[1,0;1,1;1,2;1,3];
y=[1;4;2;6];
d=zeros(4,4);
d(1,1)=1;
d(2,2)=2;
d(3,3)=3;
d(4,4)=sqrt(6);

a=(M'*M)\(M'*y);
a2=(M'*d*M)\(M'*d*y);
hold on;
plot(x,y,"*")
plot(x,a(1)+x*a(2))
plot(x,a2(1)+x*a2(2))
legend("Points","Regular","Weighted")

