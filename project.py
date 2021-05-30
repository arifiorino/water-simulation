import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 8

u = np.zeros((N+1, N))
v = np.zeros((N, N+1))
p = np.zeros((N, N))
dt = 0.01
SOLID, FLUID, EMPTY = 0, 1, 2
types = np.zeros((N, N), dtype=int)+EMPTY
for i in range(N):
  types[0][i]=SOLID
  types[i][0]=SOLID
  types[N-1][i]=SOLID
  types[i][N-1]=SOLID
for i in range(1,N-1):
  for j in range(1,N//2):
    types[i][j]=FLUID
    u[i][j]=np.random.rand()-0.5
    u[i+1][j]=np.random.rand()-0.5
    v[i][j]=np.random.rand()-0.5
    v[i][j+1]=np.random.rand()-0.5


def project():
  global u, v
  fluid_cells = []
  fluid_dict = dict()
  for i in range(N):
    for j in range(N):
      if types[i][j]==FLUID:
        fluid_cells.append((i,j))
        fluid_dict[(i,j)] = len(fluid_cells)-1
  n_fluid = len(fluid_cells)
  A = np.zeros((n_fluid, n_fluid))
  b = np.zeros((n_fluid,))
  for idx in range(n_fluid):
    i,j = fluid_cells[idx]
    if types[i-1][j]==FLUID:
      A[idx][fluid_dict[(i-1,j)]]=-1
    if types[i][j-1]==FLUID:
      A[idx][fluid_dict[(i,j-1)]]=-1
    if types[i+1][j]==FLUID:
      A[idx][fluid_dict[(i+1,j)]]=-1
    if types[i][j+1]==FLUID:
      A[idx][fluid_dict[(i,j+1)]]=-1
    non_solid = 0
    D=0
    if types[i-1][j]!=SOLID:
      D-=u[i][j]
      non_solid+=1
    if types[i][j-1]!=SOLID:
      D-=v[i][j]
      non_solid+=1
    if types[i+1][j]!=SOLID:
      D+=u[i+1][j]
      non_solid+=1
    if types[i][j+1]!=SOLID:
      D+=v[i][j+1]
      non_solid+=1
    A[idx][idx]=non_solid
    b[idx]=-(1/dt)*D
  print('eig',np.linalg.eigvals(A))
  x, _ = scipy.sparse.linalg.cg(A, b)
  for i in range(N):
    for j in range(N):
      if (i,j) in fluid_cells:
        p[i][j] = x[fluid_dict[(i,j)]]
  u2 = np.copy(u)
  for i in range(N-1):
    for j in range(N):
      #updating u[i+1][j]
      if types[i][j]!=FLUID and types[i+1][j]!=FLUID:
        continue
      rightP, leftP = 0, 0
      if types[i+1][j]==FLUID:
        rightP = x[fluid_dict[(i+1,j)]]
      if types[i+1][j]==EMPTY:
        rightP = 0
      if types[i+1][j]==SOLID:
        u2[i+1][j]=0
        continue
      if types[i][j]==FLUID:
        leftP = x[fluid_dict[(i,j)]]
      if types[i][j]==EMPTY:
        leftP = 0
      if types[i][j]==SOLID:
        u2[i+1][j]=0
        continue
      u2[i+1][j]-=dt*(rightP - leftP)
  v2 = np.copy(v)
  for i in range(N):
    for j in range(N-1):
      #updating v[i][j+1]
      if types[i][j]!=FLUID and types[i][j+1]!=FLUID:
        continue
      topP, bottomP = 0, 0
      if types[i][j+1]==FLUID:
        topP = x[fluid_dict[(i,j+1)]]
      if types[i][j+1]==EMPTY:
        topP = 0
      if types[i][j+1]==SOLID:
        v2[i][j+1]=0
        continue
      if types[i][j]==FLUID:
        bottomP = x[fluid_dict[(i,j)]]
      if types[i][j]==EMPTY:
        bottomP = 0
      if types[i][j]==SOLID:
        v2[i][j+1]=0
        continue
      v2[i][j+1]-=dt*(topP - bottomP)
  u = u2
  v = v2
  for i in range(N):
    for j in range(N):
      if types[i][j]==FLUID:
        D = u[i+1][j]-u[i][j]+v[i][j+1]-v[i][j]
        print(i,j,D)

def plot(i=0):
  plt.figure(i)
  plt.pcolormesh(p.T)
  plt.colorbar()
  plt.figure(i, figsize=(7, 7))
  plt.quiver(np.arange(0, N+0.5), np.arange(0.5, N),
             u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(0, N+0.5),
             np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()
  plt.pause(0.1)

plot(0)
project()
plot(1)

plt.show()
