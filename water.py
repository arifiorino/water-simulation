import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

# Zeroing boundary velocities
# Metaballs to render
# Smaller timestep?

N = 16

u = np.zeros((N+1, N))
v = np.zeros((N, N+1))
dt = 0.2
gy = -0.4
SOLID, FLUID, EMPTY = 0, 1, 2
types = np.zeros((N, N), dtype=int)+EMPTY
for i in range(N):
  types[0][i]=SOLID
  types[i][0]=SOLID
  types[N-1][i]=SOLID
  types[i][N-1]=SOLID
for j in range(1,N-1):
  for i in range(1,N//2):
    types[i][j]=FLUID
particles = []
for i in range(N):
  for j in range(N):
    if types[i][j]==FLUID:
      r = lambda: np.random.rand(2,)*0.5-0.25
      particles.append(np.array([i+0.25,j+0.25])+r())
      particles.append(np.array([i+0.25,j+0.75])+r())
      particles.append(np.array([i+0.75,j+0.25])+r())
      particles.append(np.array([i+0.75,j+0.75])+r())

def bilinear_interp_u(p):
  x, y = p
  x=max(x,0)
  x=min(x,N-0.0001)
  y=max(y,0.5)
  y=min(y,N-0.5001)
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
  x=max(x,0.5)
  x=min(x,N-0.5001)
  y=max(y,0)
  y=min(y,N-0.0001)
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
  k2 = np.array([bilinear_interp_u(curr + 0.5*dt*k1),
                 bilinear_interp_v(curr + 0.5*dt*k1)])
  return curr + dt * k2
def backwards_RK2(curr):
  k1 = np.array([bilinear_interp_u(curr),
                 bilinear_interp_v(curr)])
  k2 = np.array([bilinear_interp_u(curr - 0.5*dt*k1),
                 bilinear_interp_v(curr - 0.5*dt*k1)])
  return curr - dt * k2

def extrapolate():
  maxint = N*N
  d=np.zeros((N+1,N),dtype=int)+maxint
  W=[]
  for i in range(N-1):
    for j in range(N):
      if types[i][j]==FLUID or types[i+1][j]==FLUID:
        d[i+1][j]=0
  for i in range(N+1):
    for j in range(N):
      if d[i][j]==0:
        continue
      for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
        if i2 in range(N+1) and j2 in range(N):
          if d[i2][j2]==0:
            d[i][j]=1
            W.append((i,j))
            break
  t=0
  while t<len(W):
    i,j=W[t]
    new = 0
    count = 0
    for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
      if i2 in range(N+1) and j2 in range(N):
        if d[i2][j2]<d[i][j]:
          new += u[i2][j2]
          count+=1
    u[i][j]=new/count
    for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
      if i2 in range(N+1) and j2 in range(N):
        if d[i2][j2]==maxint:
          d[i2][j2]=d[i][j]+1
          W.append((i2,j2))
    t+=1
  d=np.zeros((N,N+1),dtype=int)+maxint
  W=[]
  for i in range(N):
    for j in range(N-1):
      if types[i][j]==FLUID or types[i][j+1]==FLUID:
        d[i][j+1]=0
  for i in range(N):
    for j in range(N+1):
      if d[i][j]==0:
        continue
      for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
        if i2 in range(N) and j2 in range(N+1):
          if d[i2][j2]==0:
            d[i][j]=1
            W.append((i,j))
            break
  t=0
  while t<len(W):
    i,j=W[t]
    new = 0
    count = 0
    for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
      if i2 in range(N) and j2 in range(N+1):
        if d[i2][j2]<d[i][j]:
          new += v[i2][j2]
          count+=1
    v[i][j]=new/count
    for i2,j2 in [(i-1,j),(i,j-1),(i+1,j),(i,j+1)]:
      if i2 in range(N) and j2 in range(N+1):
        if d[i2][j2]==maxint:
          d[i2][j2]=d[i][j]+1
          W.append((i2,j2))
    t+=1

def advect():
  global u, v
  u2 = np.zeros((N+1,N))
  for i in range(N+1):
    for j in range(N):
      curr = np.array([i, j+0.5])
      prev = backwards_RK2(curr)
      u2[i][j]=bilinear_interp_u(prev)
  v2 = np.zeros((N,N+1))
  for i in range(N):
    for j in range(N+1):
      curr = np.array([i+0.5, j])
      prev = backwards_RK2(curr)
      v2[i][j]=bilinear_interp_v(prev)
  u = u2
  v = v2

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
  x, _ = scipy.sparse.linalg.cg(A, b)
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

def move_particles():
  for i in range(N):
    for j in range(N):
      if types[i][j]!=SOLID:
        types[i][j]=EMPTY
  for idx in range(len(particles)):
    particles[idx]=RK2(particles[idx])
    i,j=particles[idx]
    types[int(i)][int(j)]=FLUID

def plot(i=0):
  plt.clf()
  plt.figure(i)
  plt.pcolormesh(types.T)
  plt.quiver(np.arange(0, N+0.5), np.arange(0.5, N),
             u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(0, N+0.5),
             np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")
  for x,y in particles:
    plt.plot(x, y, marker='o', markersize=5, color='black')

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()
  plt.pause(0.01)

project()
while 1:
  advect()
  project()
  extrapolate()
  move_particles()
  for i in range(N):
    for j in range(N-1):
      if types[i][j]==FLUID or types[i][j+1]==FLUID:
        v[i][j+1] += dt*gy
  plot()

plt.show()
