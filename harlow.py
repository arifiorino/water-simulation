import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 8

p = np.zeros((N, N))
u = np.random.rand(N-1, N)*0.6-0.3
v = np.random.rand(N, N-1)*0.6-0.3
#u = np.zeros((N-1,N))
#v = np.zeros((N, N-1))
BOUNDARY, FULL, SURFACE, EMPTY = 0, 1, 2, 3
types = np.zeros((N,N), dtype=int)
particles = []
dt = 0.5
dx = 1
dy = 1
atmP = 0.1 #Atmospheric pressure

gx = 0
gy = 0#-.05

def interp(f, p1, p2):
  #if f==1:
    #return interp2((p2[0]-p1[0])/2.0, (p2[1]-p1[1])/2.0)[0]
  #return interp2((p2[0]-p1[0])/2.0, (p2[1]-p1[1])/2.0)[1]
  a, b = None, None
  if p1[0]>=0 and p1[0]<len(f) and p1[1]>=0 and p1[1]<len(f[0]):
    a = f[p1[0]][p1[1]]
  if p2[0]>=0 and p2[0]<len(f) and p2[1]>=0 and p2[1]<len(f[0]):
    b = f[p2[0]][p2[1]]
  if a==None and b==None:
    return 0
  if a==None:
    return b/2.0
  if b==None:
    return a/2.0
  return (a+b)/2.0

def interp2(x, y):
  right=x+0.5
  left=x-0.5
  top=y+0.5
  bottom=y-0.5
  horizontal = int(y-0.5)+1.0
  vertical = int(x)+0.5
  TL = (vertical - left)*(top - horizontal)
  TR = (right - vertical)*(top - horizontal)
  BL = (vertical - left)*(horizontal - bottom)
  BR = (right - vertical)*(horizontal - bottom)
  i = int(vertical)
  j = int(horizontal+0.1)
  uK = TL * u[i-1][j] + TR * u[i][j] + BL * u[i-1][j-1] + BR * u[i][j-1]

  horizontal = int(y)+0.5
  vertical = int(x-0.5)+1.0
  TL = (vertical - left)*(top - horizontal)
  TR = (right - vertical)*(top - horizontal)
  BL = (vertical - left)*(horizontal - bottom)
  BR = (right - vertical)*(horizontal - bottom)
  i = int(vertical+0.1)
  j = int(horizontal)
  vK = TL * v[i-1][j] + TR * v[i][j] + BL * v[i-1][j-1] + BR * v[i][j-1]
  return (uK, vK)

def setBoundarySurface():
  global types
  #Set types
  types = np.zeros((N,N), dtype=int)+EMPTY
  for x,y in particles:
    types[int(x)][int(y)] = FULL
  for i in range(N):
    types[0][i]=BOUNDARY
    types[i][0]=BOUNDARY
    types[N-1][i]=BOUNDARY
    types[i][N-1]=BOUNDARY
  for i in range(1,N-1):
    for j in range(1,N-1):
      top = types[i][j+1] == EMPTY
      right = types[i+1][j] == EMPTY
      bottom = types[i][j-1] == EMPTY
      left = types[i-1][j] == EMPTY
      if types[i][j]==FULL and (left or right or top or bottom):
        types[i][j]=SURFACE
  #Set border vel
  for i in range(N):
    u[0][i] = 0.0
    u[N-2][i] = 0.0
    v[i][0] = 0.0
    v[i][N-2] = 0.0
  for i in range(N-1):
    u[i][0] = u[i][1]
    u[i][N-1] = u[i][N-2]
    v[0][i] = v[1][i]
    v[N-1][i] = v[N-2][i]
  #Set border pressure
  for i in range(1,N-1):
    p[0][i]=p[1][i]
    p[i][0]=p[i][1]
    p[N-1][i]=p[N-2][i]
    p[i][N-1]=p[i][N-2]
  #Set empty pressure, vel
  for i in range(N):
    for j in range(N):
      if types[i][j]==EMPTY:
        u[i][j]=0
        u[i-1][j]=0
        v[i][j]=0
        v[i][j-1]=0
        p[i][j]=atmP

def setFullPressure():
  A = np.zeros((N*N, N*N))
  b = np.zeros((N*N,))
  for i in range(N):
    for j in range(N):
      if types[i][j]!=FULL:
        continue
      UCenter = interp(u, (i,j), (i-1,j))
      URCenter = interp(u, (i,j), (i+1,j))
      ULCenter = interp(u, (i-1,j), (i-2,j))

      UTRCorner = interp(u, (i,j), (i,j+1))
      UBRCorner = interp(u, (i,j), (i,j-1))
      UTLCorner = interp(u, (i-1,j), (i-1,j+1))
      UBLCorner = interp(u, (i-1,j), (i-1,j-1))

      VCenter = interp(v, (i,j), (i,j-1))
      VTCenter = interp(v, (i,j), (i,j+1))
      VBCenter = interp(v, (i,j-1), (i,j-2))

      VTRCorner = interp(v, (i,j), (i+1,j))
      VTLCorner = interp(v, (i,j), (i-1,j))
      VBRCorner = interp(v, (i,j-1), (i+1,j-1))
      VBLCorner = interp(v, (i,j-1), (i-1,j-1))

      Q = (1/(dx**2)) * (URCenter**2 + ULCenter**2 + 2*(UCenter**2)) +\
          (1/(dy**2)) * (VTCenter**2 + VBCenter**2 + 2*(VCenter**2)) +\
          2 / (dx*dy) *\
          ((UTRCorner*VTRCorner) + (UBLCorner*VBLCorner) -\
           (UBRCorner*VBRCorner) - (UTLCorner*VTLCorner))
      D = (-u[i-1][j] + u[i][j])/dx +(-v[i][j-1] + v[i][j])/dy
      R = Q - (D/dt)
        
      A[i*N+j][(i+1)*N+j]=-1/(dx**2)
      A[i*N+j][(i-1)*N+j]=-1/(dx**2)
      A[i*N+j][i*N+j]=2/(dx**2)+2/(dy**2)
      A[i*N+j][i*N+(j+1)]=-1/(dy**2)
      A[i*N+j][i*N+(j-1)]=-1/(dy**2)
      b[i*N+j]=R
  x, info = scipy.sparse.linalg.cg(A, b)
  for i in range(N):
    for j in range(N):
      p[i][j]=x[i*N+j]

def setFullVel():
  global u, v
  u2 = np.copy(u)
  v2 = np.copy(v)
  for i in range(N-1):
    for j in range(N):
      if (types[i][j]==FULL and types[i+1][j]!=BOUNDARY) or
         (types[i+1][j]==FULL and types[i][j]!=BOUNDARY):
        UCenter = interp(u, (i,j), (i-1,j))
        URCenter = interp(u, (i,j), (i+1,j))
        UTRCorner = interp(u, (i,j), (i,j+1))
        VTRCorner = interp(v, (i,j), (i+1,j))
        UBRCorner = interp(u, (i,j), (i,j-1))
        VBRCorner = interp(v, (i,j-1), (i+1,j-1))
        PCenter = p[i][j]
        PRCenter = p[i+1][j]
        du = (1/dx)*(UCenter**2 - URCenter**2) +\
             (1/dy)*(UBRCorner*VBRCorner - UTRCorner*VTRCorner) +\
             gx +\
             (1/dx)*(PCenter - PRCenter)
        u2[i][j] += dt*du

  for i in range(N):
    for j in range(N-1):
      if (types[i][j]==FULL and types[i][j+1]!=BOUNDARY) or
         (types[i][j+1]==FULL and types[i][j]!=BOUNDARY):
        VCenter = interp(v, (i,j), (i,j-1))
        VTCenter = interp(v, (i,j), (i,j+1))
        VTRCorner = interp(v, (i,j), (i+1,j))
        UTRCorner = interp(u, (i,j), (i,j+1))
        VTLCorner = interp(v, (i,j), (i-1,j))
        UTLCorner = interp(u, (i-1,j), (i-1,j+1))
        PCenter = p[i][j]
        PTCenter = p[i][j+1]
        dv = (1/dy)*(VCenter**2 - VTCenter**2) +\
             (1/dx)*(VTLCorner*UTLCorner - VTRCorner*UTRCorner) +\
             gy +\
             (1/dy)*(PCenter - PTCenter)
        v2[i][j] += dt*dv
  u = u2
  v = v2

def setSurface():
  for i in range(N):
    for j in range(N):
      if types[i][j]!=SURFACE:
        continue
      

def moveParticles():
  for i in range(len(particles)):
    x,y = particles[i]
    uK, vK = interp2(x, y)
    particles[i] = (x+dt*uK, y+dt*vK)

def plotAll(i):
  plt.figure(i, figsize=(10, 10))
  p_types = np.zeros((N, N*2))
  k = np.amax(p)/3
  for i in range(N):
    for j in range(2*N):
      if j%2==0:
        p_types[i][j]=p[j//2][i]
      else:
        p_types[i][j]=types[(j-1)//2][i]*k
  #plt.pcolormesh(np.arange(0, N+0.1, 0.5), np.arange(N+0.1), p_types)
  plt.pcolormesh(p.T)
  plt.colorbar()
  plt.quiver(np.arange(1, N), np.arange(0.5, N), u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(1, N), np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()

def plotParticles():
  for x, y in particles:
    dx, dy = interp2(x, y)
    plt.arrow(x, y, dx, dy, head_width=0.05, head_length=0.1,
              color="black", length_includes_head=True)

for i in range(1, N-1):
  for j in range(4, 7):
    types[i][j]=EMPTY
for i in range(1, N-1):
  for j in range(1, 4):
    types[i][j]=FULL

for i in range(N):
  for j in range(N):
    if types[i][j]==FULL:
      particles.append((i+0.25, j+0.25))
      particles.append((i+0.25, j+0.75))
      particles.append((i+0.75, j+0.25))
      particles.append((i+0.75, j+0.75))

while 1:
  fixed()

  plt.clf()
  plotAll(1)
  plotParticles()
  plt.pause(0.1)

  setBoundarySurface()
  setFullPressure()
  setFullVel()
  setSurface()
  moveParticles()

plt.show()
