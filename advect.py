import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 4

u = np.random.rand(N+1, N)-0.5
v = np.random.rand(N, N+1)-0.5
dt = 0.01

def vS(i, j):
  if i<0 or i>=N or j<0 or j>=N+1:
    return 0
  return v[i][j]
def uS(i, j):
  if i<0 or i>=N+1 or j<0 or j>=N:
    return 0
  return u[i][j]
def pS(i, j):
  if i<0 or i>=N or j<0 or j>=N:
    return 0
  return p[i][j]
def typesS(i, j):
  if i<0 or i>=N or j<0 or j>=N:
    return BOUNDARY
  return types[i][j]

def bilinear_interp_u(p):
  x, y = p
  x1 = int(x)
  x2 = x1+1
  y1 = int(y-0.5)+0.5
  y2 = y1+1
  i = int(x)
  j = int(y-0.5)
  Q11 = uS(i,j)
  Q12 = uS(i,j+1)
  Q21 = uS(i+1,j)
  Q22 = uS(i+1,j+1)
  a = np.array([x2-x, x-x1])
  b = np.array([[Q11, Q12], [Q21, Q22]])
  c = np.array([y2-y, y-y1])
  return a @ b @ c.T
def bilinear_interp_v(p):
  x, y = p
  x1 = int(x-0.5)+0.5
  x2 = x1+1
  y1 = int(y)
  y2 = y1+1
  i = int(x-0.5)
  j = int(y)
  Q11 = vS(i,j)
  Q12 = vS(i,j+1)
  Q21 = vS(i+1,j)
  Q22 = vS(i+1,j+1)
  a = np.array([x2-x, x-x1])
  b = np.array([[Q11, Q12], [Q21, Q22]])
  c = np.array([y2-y, y-y1])
  return a @ b @ c.T

def advect():
  global u, v
  u2 = np.zeros((N+1, N))
  for i in range(N+1):
    for j in range(N):
      curr = np.array([i, j+0.5])
      vel = np.array([uS(i,j), (vS(i-1,j-1)+vS(i+1,j-1)+vS(i-1,j+1)+vS(i+1,j+1))/4])
      prev = curr - vel
      u2[i][j] = bilinear_interp_u(prev)
  v2 = np.zeros((N, N+1))
  for i in range(N):
    for j in range(N+1):
      curr = np.array([i+0.5, j])
      vel = np.array([(uS(i-1,j-1)+uS(i+1,j-1)+uS(i-1,j+1)+uS(i+1,j+1))/4, vS(i,j)])
      prev = curr - vel
      v2[i][j] = bilinear_interp_v(prev)
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
