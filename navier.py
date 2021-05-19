import numpy as np
import matplotlib.pyplot as plt

N = 20

p = np.zeros((N, N))
u = np.random.rand(N-1, N)-0.5
v = np.random.rand(N, N-1)-0.5
EMPTY, FULL, SURFACE, BOUNDARY = 0, 1, 2, 3
types = np.zeros((N,N), dtype=int)
particles = []
dt = 0.5
dx = 1
dy = 1
beta0 = 1.7
D_epsilon = 0.0001

gx = 0
gy = -.05

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

def moveParticles():
  for i in range(len(particles)):
    x,y = particles[i]
    uK, vK = interp2(x, y)
    particles[i] = (x+uK, y+vK)

def classify():
  types = np.zeros((N,N), dtype=int)
  for x,y in particles:
    types[int(x)][int(y)] = FULL
  for i in range(1,N-1):
    for j in range(1,N-1):
      top = types[i][j+1] == EMPTY
      right = types[i+1][j] == EMPTY
      bottom = types[i][j-1] == EMPTY
      left = types[i-1][j] == EMPTY
      if types[i][j]==FULL and (left or right or top or bottom):
        types[i][j]=SURFACE
  for i in range(N):#Border is boundary
    types[0][i]=BOUNDARY
    types[i][0]=BOUNDARY
    types[N-1][i]=BOUNDARY
    types[i][N-1]=BOUNDARY

def setBoundary():
  for i in range(N):
    u[0][i] = 0.0
    u[N-2][i] = 0.0
    v[i][0] = 0.0
    v[i][N-2] = 0.0
  for i in range(N-1):
    u[i][0] = 0.0
    u[i][N-1] = 0.0
    v[0][i] = 0.0
    v[N-1][i] = 0.0
  
def setSurface():
  for i in range(N):
    for j in range(N):
      if types[i][j]!=SURFACE:
        continue
      into = 0.0
      count = 4
      if types[i-1][j]==FULL:
        into += u[i-1][j]
        count-=1
      if types[i][j-1]==FULL:
        into += v[i][j-1]
        count-=1
      if types[i+1][j]==FULL:
        into -= u[i][j]
        count-=1
      if types[i][j+1]==FULL:
        into -= v[i][j]
        count-=1
      if types[i-1][j]==EMPTY:
        u[i-1][j]=-into/count
      if types[i][j-1]==EMPTY:
        v[i][j-1]=-into/count
      if types[i+1][j]==EMPTY:
        u[i][j]=into/count
      if types[i][j+1]==EMPTY:
        v[i][j]=into/count
 
def navier():
  for i in range(N):
    for j in range(N):
      if types[i][j]!=FULL:
        continue
      UCenter = interp(u, (i,j), (i-1,j))
      URCenter = interp(u, (i,j), (i+1,j))
      UTRCorner = interp(u, (i,j), (i,j+1))
      VTRCorner = interp(v, (i,j), (i,j+1))
      UBRCorner = interp(u, (i,j), (i,j-1))
      VBRCorner = interp(v, (i,j), (i,j-1))
      PCenter = p[i][j]
      PRCenter = p[i+1][j]
      du = (1/dx)*(UCenter**2 - URCenter**2) +\
           (1/dy)*(UBRCorner*VBRCorner - UTRCorner*VTRCorner) +\
           gx +\
           (1/dx)*(PCenter - PRCenter)
      u[i][j] += dt*du

      VCenter = interp(v, (i,j), (i,j-1))
      VTCenter = interp(v, (i,j), (i,j+1))
      VTRCorner = interp(v, (i,j), (i+1,j))
      UTRCorner = interp(u, (i,j), (i+1,j))
      VTLCorner = interp(v, (i,j), (i-1,j))
      UTLCorner = interp(u, (i,j), (i-1,j))
      PCenter = p[i][j]
      PTCenter = p[i][j+1]
      dv = (1/dy)*(VCenter**2 - VTCenter**2) +\
           (1/dx)*(VTLCorner*UTLCorner - VTRCorner*UTRCorner) +\
           gy +\
           (1/dy)*(PCenter - PTCenter)
      v[i][j] += dt*dv
      

def pressureSolve():
  worstD = D_epsilon+1
  beta = beta0 / 2.0 * dt *(1.0/(dx*dx) + 1.0/(dy*dy))
  count = 0
  while worstD > D_epsilon:
    worstD = 0
    for i in range(N):
      for j in range(N):
        if types[i][j] != FULL:
          continue
        D = u[i-1][j] - u[i][j] + v[i][j-1] - v[i][j]
        if abs(D)>worstD:
          worstD = abs(D)
        dp = beta * D
        u[i][j]+=(dt/dx)*dp
        u[i-1][j]-=(dt/dx)*dp
        v[i][j]+=(dt/dy)*dp
        v[i-1][j]-=(dt/dy)*dp
        p[i][j]+=dp
    count+=1
  #print(count, worstD)

def plotAll(i,center=False):
  plt.figure(i, figsize=(7, 7))
  plt.quiver(np.arange(1, N), np.arange(0.5, N), u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1, color="blue")
  plt.quiver(np.arange(0.5, N), np.arange(1, N), np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1, color="red")

  if center:
    cx = np.zeros((N-2,N-2))
    cy = np.zeros((N-2,N-2))
    for i in range(1,N-1):
      for j in range(1,N-1):
        cx[i-1][j-1]=(u[i][j]+u[i-1][j])/2.0
        cy[i-1][j-1]=(v[i][j-1]+v[i][j])/2.0
    plt.quiver(np.arange(1.5, N-1), np.arange(1.5, N-1), cx.T, cy.T, angles='xy', scale_units='xy', scale=1, color="pink")

  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()

def plotParticles():
  for x, y in particles:
    dx, dy = interp2(x, y)
    plt.arrow(x, y, dx, dy, head_width=0.05, head_length=0.1,
              color="green", length_includes_head=True)


for i in range(1, N-1):
  for j in range(1, N-1):
    types[i][j]=FULL

plt.show()
