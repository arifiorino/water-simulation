import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 16

u = np.zeros((N+1, N))
v = np.zeros((N, N+1))
dt = 0.2
gy = -0.4
SOLID, FLUID, EMPTY = 0, 1, 2
types = np.zeros((N, N), dtype=int)+EMPTY
phi = np.zeros((N, N))
for i in range(N):
  types[0][i]=SOLID
  types[i][0]=SOLID
  types[N-1][i]=SOLID
  types[i][N-1]=SOLID
for i in range(1,N//2):
  for j in range(1,N-1):
    types[i][j]=FLUID
for i in range(N):
  for j in range(N):
    if types[i][j]==FLUID:
      phi[i][j]=-1
    else:
      phi[i][j]=1
for i in range(0,N//2):
  phi[i][0]=-1
  phi[i][N-1]=-1
  phi[0][i]=-1
  phi[0][N-1-i]=-1

def calc_phi():
  global phi
  phi2 = np.zeros((N, N)) 
  findZero = lambda y1, y2: -y1 / (y2-y1)
  closest_point = [[None for i in range(N)] for j in range(N)]
  known = [[False for i in range(N)] for j in range(N)]
  for i in range(N):
    for j in range(N):
      if types[i][j]==SOLID:
        continue
      distances=[]
      candidates=[]
      if i<N-1 and (phi[i][j]<0) != (phi[i+1][j]<0) and types[i+1][j]!=SOLID:
        d=findZero(phi[i][j],phi[i+1][j])
        distances.append(d)
        candidates.append((i+0.5+d,j+0.5))
      if i>0 and (phi[i][j]<0) != (phi[i-1][j]<0) and types[i-1][j]!=SOLID:
        d=findZero(phi[i][j],phi[i-1][j])
        distances.append(d)
        candidates.append((i+0.5-d,j+0.5))
      if j<N-1 and (phi[i][j]<0) != (phi[i][j+1]<0) and types[i][j+1]!=SOLID:
        d=findZero(phi[i][j],phi[i][j+1])
        distances.append(d)
        candidates.append((i+0.5,j+0.5+d))
      if j>0 and (phi[i][j]<0) != (phi[i][j-1]<0) and types[i][j-1]!=SOLID:
        d=findZero(phi[i][j],phi[i][j-1])
        distances.append(d)
        candidates.append((i+0.5,j+0.5-d))
      if len(candidates)>0:
        known[i][j]=True
        closest_point[i][j]=candidates[np.argmin(distances)]
        phi2[i][j]=min(distances)
  def loop_for(i,j):
    if known[i][j]:
       return
    distances = []
    candidates=[]
    if i<N-1 and known[i+1][j]:
      d = (i+0.5 - closest_point[i+1][j][0])**2 + (j+0.5 - closest_point[i+1][j][1])**2
      candidates.append(closest_point[i+1][j])
      distances.append(d)
    if i>0 and known[i-1][j]:
      d = (i+0.5 - closest_point[i-1][j][0])**2 + (j+0.5 - closest_point[i-1][j][1])**2
      candidates.append(closest_point[i-1][j])
      distances.append(d)
    if j<N-1 and known[i][j+1]:
      d = (i+0.5 - closest_point[i][j+1][0])**2 + (j+0.5 - closest_point[i][j+1][1])**2
      candidates.append(closest_point[i][j+1])
      distances.append(d)
    if j>0 and known[i][j-1]:
      d = (i+0.5 - closest_point[i][j-1][0])**2 + (j+0.5 - closest_point[i][j-1][1])**2
      candidates.append(closest_point[i][j-1])
      distances.append(d)
    if len(distances)==0:
      return
    if closest_point[i][j]!=None:
      candidates.append(closest_point[i][j])
      distances.append(phi2[i][j]**2)

    closest_point[i][j]=candidates[np.argmin(distances)]
    phi2[i][j]=np.sqrt(min(distances))
    known[i][j]=True
    if i<N-1 and known[i+1][j] and phi2[i+1][j]>phi2[i][j]:
      known[i+1][j]=False
    if i>0 and known[i-1][j] and phi2[i-1][j]>phi2[i][j]:
      known[i-1][j]=False
    if j<N-1 and known[i][j+1] and phi2[i][j+1]>phi2[i][j]:
      known[i][j+1]=False
    if j>0 and known[i][j-1] and phi2[i][j-1]>phi2[i][j]:
      known[i][j-1]=False
  for _ in range(2):
    for j in range(N):
      for i in range(N):
        loop_for(i,j)
    for j in range(N):
      for i in range(N-1,-1,-1):
        loop_for(i,j)
    for i in range(N):
      for j in range(N):
        loop_for(i,j)
    for i in range(N):
      for j in range(N-1,-1,-1):
        loop_for(i,j)
  for i in range(N):
    for j in range(N):
      if phi[i][j]<0:
        phi2[i][j]*=-1
  phi = phi2
  for i in range(N):
    for j in range(N):
      if types[i][j]!=SOLID:
        if phi[i][j]<0:
          types[i][j]=FLUID
        else:
          types[i][j]=EMPTY

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
def bilinear_interp_phi(p):
  x, y = p
  x=max(x,0.5)
  x=min(x,N-0.5001)
  y=max(y,0.5)
  y=min(y,N-0.5001)
  x1 = int(x-0.5)+0.5
  x2 = x1+1
  y1 = int(y-0.5)+0.5
  y2 = y1+1
  i = int(x-0.5)
  j = int(y-0.5)
  Q11 = phi[i][j]
  Q12 = phi[i][j+1]
  Q21 = phi[i+1][j]
  Q22 = phi[i+1][j+1]
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

def advect_phi():
  global phi
  phi2 = np.zeros((N,N))
  for i in range(N):
    for j in range(N):
      curr = np.array([i+0.5, j+0.5])
      prev = backwards_RK2(curr)
      phi2[i][j]=bilinear_interp_phi(prev)
  phi=phi2

def advect_uv():
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

def plot(i=0):
  plt.clf()
  plt.figure(i)
  colors=np.zeros((N,N))
  for i in range(N):
    for j in range(N):
      if phi[i][j]<0:
        colors[i][j]=phi[i][j]
  plt.pcolormesh(colors.T)
  plt.colorbar()
  plt.quiver(np.arange(0, N+0.5), np.arange(0.5, N),
             u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(0, N+0.5),
             np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()
  plt.pause(0.01)

project()
calc_phi()
while 1:
  calc_phi()
  for i in range(N):
    for j in range(N-1):
      if types[i][j]==EMPTY and types[i][j+1]==EMPTY:
        v[i][j+1]=0
  for i in range(N-1):
    for j in range(N):
      if types[i][j]==EMPTY and types[i+1][j]==EMPTY:
        u[i][j+1]=0
  for i in range(N):
    for j in range(N-1):
      if types[i][j]==FLUID or types[i][j+1]==FLUID:
        v[i][j+1] += dt*gy
  advect_uv()
  advect_phi()
  project()
  plot()

