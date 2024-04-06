import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate




f= lambda x: np.sin(9*x)


n=40
x=np.linspace(-1,1,n)
y=f(x)
points=np.linspace(-1, 1,10001)
actual=f(points)
cs=interpolate.CubicSpline(x,y)
spline=cs(points)
error=abs(spline-actual)

plt.plot(x,y,"o")
plt.plot(points,spline)
plt.plot(points,actual)
plt.show()
plt.plot(points,np.log10(error))
plt.show

