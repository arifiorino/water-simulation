import numpy as np
import matplotlib.pyplot as plt

count=50
h = 0.1 #dt

actual = lambda t : np.exp((t-np.sin(t)*np.cos(t))/2)
f = lambda t, y: np.sin(t)*np.sin(t)*y #dy/dt
y0 = 1

ts = np.zeros((count,))
ys = np.zeros((count,))

t=0
y=y0
for n in range(count):
  ts[n]=t
  ys[n]=y
  k1 = f(t,y)
  y+=h*k1
  t+=h
print('euler error',y,actual(t))
plt.plot(ts, ys, color='green')

t=0
y=y0
for n in range(count):
  ts[n]=t
  ys[n]=y
  k1 = f(t,y)
  k2 = f(t+h/2, y+h*(k1/2))
  y+=h*k2
  t+=h
print('RK2 error',y,actual(t))
plt.plot(ts, ys, color='orange')

t=0
y=y0
for n in range(count):
  ts[n]=t
  ys[n]=y
  k1 = f(t,y)
  k2 = f(t+h/2, y+h*(k1/2))
  k3 = f(t+h/2, y+h*(k2/2))
  k4 = f(t+h, y+h*k3)
  y+=(1/6)*h*(k1+2*k2+2*k3+k4)
  t+=h
print('RK4 error',y,actual(t))
plt.plot(ts, ys, color='red')

print('blue actual, green euler, orange RK2, red RK4')
plt.plot(ts, actual(ts), color='blue')
plt.show()
