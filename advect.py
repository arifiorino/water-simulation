import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 8

#u = np.random.rand(N+1, N)-0.5
#v = np.random.rand(N, N+1)-0.5
u = np.zeros((N+1, N))
v = np.zeros((N, N+1))
for i in range(N):
  for j in range(N//2+1):
    u[j][i]=-0.5
    u[N-j][i]=0.5
    v[i][j]=-0.5
    v[i][N-j]=0.5

dt = 0.1

def bilinear_interp_u(p):
  x, y = p
  #clamp
  x=max(x,0)
  x=min(x,N-0.0001)
  y=max(y,0.5)
  y=min(y,N-0.5001)
  #from wiki
  x1 = int(x)
  x2 = x1+1
  y1 = int(y-0.5)+0.5
  y2 = y1+1
  i = int(x)
  j = int(y-0.5)
  Q11 = u[i][j]
  Q12 = u[i][j+1]
  Q21 = u[i+1][j]
  Q22 = u[i+1][j+1]
  a = np.array([x2-x, x-x1])
  b = np.array([[Q11, Q12], [Q21, Q22]])
  c = np.array([y2-y, y-y1])
  return a @ b @ c.T
def bilinear_interp_v(p):
  x, y = p
  #clamp
  x=max(x,0.5)
  x=min(x,N-0.5001)
  y=max(y,0)
  y=min(y,N-0.0001)
  #from wiki
  x1 = int(x-0.5)+0.5
  x2 = x1+1
  y1 = int(y)
  y2 = y1+1
  i = int(x-0.5)
  j = int(y)
  Q11 = v[i][j]
  Q12 = v[i][j+1]
  Q21 = v[i+1][j]
  Q22 = v[i+1][j+1]
  a = np.array([x2-x, x-x1])
  b = np.array([[Q11, Q12], [Q21, Q22]])
  c = np.array([y2-y, y-y1])
  return a @ b @ c.T

def RK2(curr):
  k1 = np.array([bilinear_interp_u(curr),
                 bilinear_interp_v(curr)])
  k2 = np.array([bilinear_interp_u(curr - 0.5*dt*k1),
                 bilinear_interp_v(curr - 0.5*dt*k1)])
  return - dt * k2

def advect():
  global u, v
  u2 = np.zeros((N+1,N))
  for i in range(N+1):
    for j in range(N):
      curr = np.array([i, j+0.5])
      prev = curr + RK2(curr)
      u2[i][j]=bilinear_interp_u(prev)
      
  v2 = np.zeros((N,N+1))
  for i in range(N):
    for j in range(N+1):
      curr = np.array([i+0.5, j])
      prev = curr + RK2(curr)
      v2[i][j]=bilinear_interp_v(prev)
  u = u2
  v = v2

def plot(i=0):
  plt.clf()
  plt.figure(i, figsize=(7, 7))
  plt.quiver(np.arange(0, N+0.5), np.arange(0.5, N),
             u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(0, N+0.5),
             np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()
  plt.pause(0.1)


while 1:
  plot()
  advect()

plt.show()
