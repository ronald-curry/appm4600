import numpy as np; # import numpy library and call it np
from matplotlib import pyplot as plt; #import pyplot module and call it plt

# question 1
x = np.linspace(1.920,2.080,160);
def coeff(x):
    return x**9-18*x**8+144*x**7-672*x**6+2016*x**5-4032*x**4+5376*x**3-4608*x**2+2304*x-512;
plt.plot(x,coeff(x)) #creates the plot figure
plt.title("1) Coefficients"); #figure title

plt.show(); #necessary when executing from command line to generate the pop-up figure
def factored(x):
    return (x-2)**9;
plt.plot(x,factored(x))
plt.title("1) Factored");
plt.show();



#Question 5
x_1=.3;
x_2=1234567;
n=np.linspace(0,16,17);
delta_range=10**(-n);
def orig5 (x,delta):
    return np.cos(x+delta)-np.cos(x);
def fixed5(x,delta):
    return -2*np.sin((2*x+delta)/2)*np.sin(delta/2);
y_1=np.zeros(17);
y_2=np.zeros(17);
for j in range(17):
    y_1[j]=orig5(x_1,delta_range[j])-fixed5(x_1,delta_range[j])
plt.plot(delta_range,y_1)
plt.title("x=0.3")
plt.xscale("log")
plt.xlabel("Delta")
plt.ylabel("Orginal-manipulated expression")
plt.show()
for j in range(17):
    y_2[j]=orig5(x_2,delta_range[j])-fixed5(x_2,delta_range[j])
plt.plot(delta_range,y_2)
plt.title("x=1234567")
plt.xscale("log")
plt.xlabel("Delta")
plt.ylabel("Orginal-manipulated expression")
plt.show()