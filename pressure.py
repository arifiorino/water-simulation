import numpy as np
import matplotlib.pyplot as plt

N = 20

p = np.zeros((N, N))
u = np.random.rand(N-1, N)-0.5
v = np.random.rand(N, N-1)-0.5
dt = 0.5
dx = 1
dy = 1
beta0 = 1.7
D_epsilon = 0.0001

def pressureSolve():
  worstD = D_epsilon+1
  beta = beta0 / 2.0 * dt *(1.0/(dx*dx) + 1.0/(dy*dy))
  count = 0
  while worstD > D_epsilon:
    worstD = 0
    for i in range(1, N-1):
      for j in range(1, N-1):
        D = u[i-1][j] - u[i][j] + v[i][j-1] - v[i][j]
        if abs(D)>worstD:
          worstD = abs(D)
        dp = beta * D
        u[i][j]+=(dt/dx)*dp
        u[i-1][j]-=(dt/dx)*dp
        v[i][j]+=(dt/dy)*dp
        v[i-1][j]-=(dt/dy)*dp
        p[i][j]+=dp
    #print(worstD)
    count+=1
  print(count, worstD)

def plotVel(i):
  plt.figure(i, figsize=(7, 7))
  plt.quiver(np.arange(1, N), np.arange(0.5, N), u.T, np.zeros(u.T.shape), angles='xy', scale_units='xy', scale=1)
  plt.quiver(np.arange(0.5, N), np.arange(1, N), np.zeros(v.T.shape), v.T, angles='xy', scale_units='xy', scale=1)

  cx = np.zeros((N-2,N-2))
  cy = np.zeros((N-2,N-2))
  for i in range(1,N-1):
    for j in range(1,N-1):
      cx[i-1][j-1]=(u[i][j]+u[i-1][j])/2.0
      cy[i-1][j-1]=(v[i][j-1]+v[i][j])/2.0
  plt.quiver(np.arange(1.5, N-1), np.arange(1.5, N-1), cx.T, cy.T, angles='xy', scale_units='xy', scale=1)
  
  plt.xticks(np.arange(0, N+1))
  plt.yticks(np.arange(0, N+1))
  plt.grid()

plotVel(1)
pressureSolve()
plotVel(2)

plt.show()
