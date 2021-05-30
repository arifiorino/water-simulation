import scipy.sparse.linalg
import numpy as np
import matplotlib.pyplot as plt

N = 8

p = np.zeros((N, N))
#u = np.random.rand(N+1, N)*0.2-0.1
#v = np.random.rand(N, N+1)*0.2-0.1
u = np.zeros((N+1,N))
v = np.zeros((N, N+1))
BOUNDARY, FULL, SURFACE, EMPTY = 0, 1, 2, 3
types = np.zeros((N,N), dtype=int)
particles = []
dt = 0.2
dx = 1
dy = 1
atmP = 0 #Atmospheric pressure
waterP = 0.2

gx = 0
gy = -.04

def interp(f, p1, p2):
  return (f[p1[0]][p1[1]]+f[p2[0]][p2[1]])/2.0

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
  j = int(horizontal-0.1)
  uK = TL * u[i][j+1] + TR * u[i+1][j+1] + BL * u[i][j] + BR * u[i+1][j]

  horizontal = int(y)+0.5
  vertical = int(x-0.5)+1.0
  TL = (vertical - left)*(top - horizontal)
  TR = (right - vertical)*(top - horizontal)
  BL = (vertical - left)*(horizontal - bottom)
  BR = (right - vertical)*(horizontal - bottom)
  i = int(vertical-0.1)
  j = int(horizontal)
  vK = TL * v[i][j+1] + TR * v[i+1][j+1] + BL * v[i][j] + BR * v[i+1][j]
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
  #Set border vel - normal
  for i in range(1,N-1):
    u[0][i] = -u[2][i]
    u[N][i] = -u[N-2][i]
    v[i][0] = -v[i][2]
    v[i][N] = -v[i][N-2]
  #Set border vel - wall
  for i in range(1,N-1):
    u[1][i] = 0
    u[N-1][i] = 0
    v[i][1] = 0
    v[i][N-1] = 0
  #Set border vel - tangential
  for i in range(1,N):
    u[i][0] = u[i][1]
    u[i][N-1] = u[i][N-2]
    v[0][i] = v[1][i]
    v[N-1][i] = v[N-2][i]
  u[0][0]=0
  v[0][0]=0
  u[N][0]=0
  v[0][N]=0
  u[0][N-1]=0
  v[N-1][0]=0
  u[N][N-1]=0
  v[N-1][N]=0

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
        u[i+1][j]=0
        u[i][j]=0
        v[i][j+1]=0
        v[i][j]=0
        p[i][j]=atmP

def setSurface():
  global u, v
  u2 = np.copy(u)
  v2 = np.copy(v)
  for i in range(N):
    for j in range(N):
      if types[i][j]==SURFACE:
        #     EMPTY
        #     L B T R
        m = {(1,0,0,0): (-v[i][j]-v[i][j+1]+u[i+1][j],None,None,None),
             (0,1,0,0): (None,-u[i][j]+v[i][j+1]+u[i+1][j],None,None),
             (0,0,1,0): (None,None,u[i][j]+v[i][j]-u[i+1][j],None),
             (0,0,0,1): (None,None,None,u[i][j]+v[i][j]-v[i][j+1]),
             
             (1,1,0,0): (u[i+1][j],v[i][j+1],None,None,None),
             (1,0,1,0): (u[i+1][j],None,v[i][j],None,None),
             (1,0,0,1): (None,None,None,None),
             (0,1,1,0): (None,None,None,None),
             (0,1,0,1): (None,v[i][j+1],None,u[i][j]),
             (0,0,1,1): (None,None,v[i][j],u[i][j]),
             
             (1,1,1,0): (u[i+1][j],None,None,None),
             (1,1,0,1): (None,v[i][j+1],None,None),
             (1,0,1,1): (None,None,v[i][j],None),
             (0,1,1,1): (None,None,None,u[i][j]),
             
             (1,1,1,1): (gx,gy,gy,gx)}
        diff = m[(int(types[i-1][j]==EMPTY), int(types[i][j-1]==EMPTY),
                  int(types[i][j+1]==EMPTY), int(types[i+1][j]==EMPTY))]
        if diff[0]:
          u2[i][j]=diff[0]
        if diff[1]:
          v2[i][j]=diff[1]
        if diff[2]:
          v2[i][j+1]=diff[2]
        if diff[3]:
          u2[i+1][j]=diff[3]

        p[i][j]=atmP
  u = u2
  v = v2

def pressureStep():
  A = np.zeros((N*N, N*N))
  b = np.zeros((N*N,))
  for i in range(N):
    for j in range(N):
      if types[i][j]!=FULL:
        continue
      UCenter = interp(u, (i,j), (i+1,j))
      ULCenter = interp(u, (i,j), (i-1,j))
      URCenter = interp(u, (i+1,j), (i+2,j))

      UTRCorner = interp(u, (i+1,j), (i+1,j+1))
      UBRCorner = interp(u, (i+1,j), (i+1,j-1))
      UTLCorner = interp(u, (i,j), (i,j+1))
      UBLCorner = interp(u, (i,j), (i,j-1))

      VCenter = interp(v, (i,j), (i,j+1))
      VBCenter = interp(v, (i,j), (i,j-1))
      VTCenter = interp(v, (i,j+1), (i,j+2))

      VTRCorner = interp(v, (i,j+1), (i+1,j+1))
      VTLCorner = interp(v, (i,j+1), (i-1,j+1))
      VBRCorner = interp(v, (i,j), (i+1,j))
      VBLCorner = interp(v, (i,j), (i-1,j))

      Q = (1/(dx**2)) * (URCenter**2 + ULCenter**2 + 2*(UCenter**2)) +\
          (1/(dy**2)) * (VTCenter**2 + VBCenter**2 + 2*(VCenter**2)) +\
          2 / (dx*dy) *\
          ((UTRCorner*VTRCorner) + (UBLCorner*VBLCorner) -\
           (UBRCorner*VBRCorner) - (UTLCorner*VTLCorner))
      D = (-u[i][j] + u[i+1][j])/dx +(-v[i][j] + v[i][j+1])/dy
      R = Q - (D/dt)
        
      b[i*N+j]=R
      if types[i+1][j]==FULL:
        A[i*N+j][(i+1)*N+j]=-1/(dx**2)
      if types[i-1][j]==FULL:
        A[i*N+j][(i-1)*N+j]=-1/(dx**2)
      A[i*N+j][i*N+j]=2/(dx**2)+2/(dy**2)
      if types[i][j+1]==FULL:
        A[i*N+j][i*N+(j+1)]=-1/(dy**2)
      if types[i][j-1]==FULL:
        A[i*N+j][i*N+(j-1)]=-1/(dy**2)
  #print(np.linalg.eigvals(A))
  #x, _, _, _ = np.linalg.lstsq(A, b)
  x, info = scipy.sparse.linalg.cg(A, b)
  for i in range(N):
    for j in range(N):
      if types[i][j]==FULL:
        p[i][j]=x[i*N+j]

def velStep():
  global u, v
  u2 = np.copy(u)
  v2 = np.copy(v)
  for i in range(N-1):
    for j in range(N):
      if (types[i][j]==FULL and types[i+1][j]!=BOUNDARY) or\
         (types[i+1][j]==FULL and types[i][j]!=BOUNDARY) or\
         (types[i][j]==SURFACE and types[i+1][j]==SURFACE) or\
         (types[i+1][j]==SURFACE and types[i][j]==SURFACE):
        UCenter = interp(u, (i,j), (i+1,j))
        URCenter = interp(u, (i+1,j), (i+2,j))
        UTRCorner = interp(u, (i+1,j), (i+1,j+1))
        VTRCorner = interp(v, (i,j+1), (i+1,j+1))
        UBRCorner = interp(u, (i,j), (i,j+1))
        VBRCorner = interp(v, (i,j), (i+1,j))
        PCenter = p[i][j]
        PRCenter = p[i+1][j]
        du = (1/dx)*(UCenter**2 - URCenter**2) +\
             (1/dy)*(UBRCorner*VBRCorner - UTRCorner*VTRCorner) +\
             gx +\
             (1/dx)*(PCenter - PRCenter)
        u2[i+1][j] += dt*du

  for i in range(N):
    for j in range(N-1):
      if (types[i][j]==FULL and types[i][j+1]!=BOUNDARY) or\
         (types[i][j+1]==FULL and types[i][j]!=BOUNDARY) or\
         (types[i][j]==SURFACE and types[i][j+1]==SURFACE) or\
         (types[i][j+1]==SURFACE and types[i][j]==SURFACE):

        VCenter = interp(v, (i,j), (i,j+1))
        VTCenter = interp(v, (i,j+1), (i,j+2))
        VTRCorner = interp(v, (i,j+1), (i+1,j+1))
        UTRCorner = interp(u, (i+1,j), (i+1,j+1))
        VTLCorner = interp(v, (i,j+1), (i-1,j+1))
        UTLCorner = interp(u, (i,j), (i,j+1))
        PCenter = p[i][j]
        PTCenter = p[i][j+1]
        dv = (1/dy)*(VCenter**2 - VTCenter**2) +\
             (1/dx)*(VTLCorner*UTLCorner - VTRCorner*UTRCorner) +\
             gy +\
             (1/dy)*(PCenter - PTCenter)
        v2[i][j+1] += dt*dv
  u = u2
  v = v2


def moveParticles():
  for i in range(len(particles)):
    x,y = particles[i]
    uK, vK = interp2(x, y)
    particles[i] = (x+dt*uK, y+dt*vK)

def plotAll(i):
  plt.figure(i, figsize=(7, 7))
  plt.pcolormesh(p.T)
  plt.colorbar()
  plt.quiver(np.arange(0, N+0.5), np.arange(0.5, N), u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="red")
  plt.quiver(np.arange(0.5, N), np.arange(0, N+0.5), np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  for x, y in particles:
    dx, dy = interp2(x, y)
    plt.arrow(x, y, dx, dy, head_width=0.05, head_length=0.1,
              color="black", length_includes_head=True)

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()


for i in range(1, N-1):
  for j in range(1, N//2):
    types[i][j]=FULL

for i in range(N):
  for j in range(N):
    if types[i][j]==FULL:
      particles.append((i+0.25, j+0.25))
      particles.append((i+0.25, j+0.75))
      particles.append((i+0.75, j+0.25))
      particles.append((i+0.75, j+0.75))
      p[i][j]=waterP

for i in range(N-1):
  for j in range(N-1):
    if types[i][j]==FULL and types[i+1][j]==FULL:
      u[i+1][j]=np.random.rand()-0.5
    if types[i][j]==FULL and types[i][j+1]==FULL:
      v[i][j+1]=np.random.rand()-0.5

while 1:
  setBoundarySurface()
  setSurface()
  pressureStep()
  velStep()

  setBoundarySurface()
  setSurface()

  moveParticles()

  plt.clf()
  plotAll(1)
  plt.pause(0.1)
  #input('update?')

plt.show()
