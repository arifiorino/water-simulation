import numpy as np
import matplotlib.pyplot as plt

N = 8

p = np.zeros((N, N))
u = np.random.rand(N-1, N)-0.5
v = np.random.rand(N, N-1)-0.5
#u = np.zeros((N-1,N))
#v = np.zeros((N, N-1))
BOUNDARY, FULL, SURFACE, EMPTY = 0, 1, 2, 3
types = np.zeros((N,N), dtype=int)
particles = []
dt = 0.3
dx = 1
dy = 1
beta0 = 1.7
D_epsilon = 0.0001

gx = 0
gy = -.05

def interp(f, p1, p2):
  a, b = None, None
  if p1[0]>=0 and p1[0]<len(f) and p1[1]>=0 and p1[1]<len(f[0]):
    a = f[p1[0]][p1[1]]
  if p2[0]>=0 and p2[0]<len(f) and p2[1]>=0 and p2[1]<len(f[0]):
    b = f[p2[0]][p2[1]]
  if a==None and b==None:
    return 0
  if a==None:
    return b
  if b==None:
    return a
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

def advectin():
  global u, v
  u2 = np.copy(u)
  v2 = np.copy(v)
  for i in range(1, N-1):
    for j in range(1, N-1):
      if types[i][j]!=FULL:
        continue
      if types[i-1][j]!=BOUNDARY:
        x, y = i, j+0.5
        for k in range(5):
          tmp = interp2(x,y)
          x-=(dt/5)*tmp[0]
          y-=(dt/5)*tmp[1]
        tmp = interp2(x,y)
        u2[i-1][j]=tmp[0]

      if types[i][j-1]!=BOUNDARY:
        x, y = i+0.5, j
        for k in range(5):
          tmp = interp2(x,y)
          x-=(dt/5)*tmp[0]
          y-=(dt/5)*tmp[1]
        tmp = interp2(x,y)
        v2[i][j-1]=tmp[0]
  u = u2
  v = v2
        

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
  global u, v
  u2 = np.copy(u)
  v2 = np.copy(v)
  for i in range(N):
    for j in range(N):
      if types[i][j]!=FULL:
        continue
      #print(i,j)
      if types[i+1][j]!=BOUNDARY:
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

      if types[i][j+1]!=BOUNDARY:
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
      

def pressureSolve():
  worstD = D_epsilon+1
  beta = 1.0 / (2.0 * dt *(1.0/(dx*dx) + 1.0/(dy*dy)))
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
        x=4
        if types[i+1][j]==BOUNDARY:
          x-=1
        if types[i-1][j]==BOUNDARY:
          x-=1
        if types[i][j+1]==BOUNDARY:
          x-=1
        if types[i][j-1]==BOUNDARY:
          x-=1

        if types[i+1][j]!=BOUNDARY:
          u[i][j]+=(dt/dx)*dp*(4/x)
        if types[i-1][j]!=BOUNDARY:
          u[i-1][j]-=(dt/dx)*dp*(4/x)
        if types[i][j+1]!=BOUNDARY:
          v[i][j]+=(dt/dy)*dp*(4/x)
        if types[i][j-1]!=BOUNDARY:
          v[i][j-1]-=(dt/dy)*dp*(4/x)
        p[i][j]+=dp
    count+=1
    #print(count)
  #print(count, worstD)

def plotAll(i,center=False):
  plt.figure(i, figsize=(7, 7))
  plt.pcolormesh(p)
  plt.colorbar()
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

setBoundary()
pressureSolve()
plotAll(1)
while 1:
  plt.clf()
  plotAll(1)
  plt.pause(0.1)
  navier()
  #setBoundary()
  pressureSolve()
'''
setBoundary()
pressureSolve()
setBoundary()
plotAll(1)
for i in range(40):
  particles.append(np.random.rand(2,)*N/2+N/4)
for i in range(40):
  plotParticles()
  moveParticles()
'''

plt.show()
