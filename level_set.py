import matplotlib.pyplot as plt
import numpy as np

N = 32

phi = np.zeros((N, N))
surfaceF = lambda x: (np.sin(2*np.pi*x/(N/2)))*N/8+N/2
surface = [(x, surfaceF(x)) for x in np.arange(0,N,0.1)]
for i in range(N):
  for j in range(N):
    p=(i+0.5, j+0.5)
    phi[i][j]=np.sqrt(min([(p[0]-surface[k][0])**2 + (p[1]-surface[k][1])**2 for k in range(len(surface))]))
    if surfaceF(i+0.5) > j+0.5:
      phi[i][j]*=-1
findZero = lambda y1, y2: -y1 / (y2-y1)

surface = []
for i in range(N):
  for j in range(N):
    if i<N-1 and (phi[i][j]<0) != (phi[i+1][j]<0):
      surface.append((i+0.5+findZero(phi[i][j],phi[i+1][j]), j+0.5))
    if j<N-1 and (phi[i][j]<0) != (phi[i][j+1]<0):
      surface.append((i+0.5, j+0.5+findZero(phi[i][j],phi[i][j+1])))

plt.figure(0)
plt.pcolormesh(phi.T)
plt.colorbar()
plt.plot(np.arange(0,N,0.1),surfaceF(np.array(np.arange(0,N,0.1))))

for i in range(N):
  for j in range(N):
    p=(i+0.5, j+0.5)
    v=np.sqrt(min([(p[0]-surface[k][0])**2 + (p[1]-surface[k][1])**2 for k in range(len(surface))]))
    if phi[i][j]<0:
      v*=-1
    phi[i][j]=v
plt.figure(1)
plt.pcolormesh(phi.T)
plt.colorbar()
plt.plot([p[0] for p in surface], [p[1] for p in surface], 'o')
plt.plot(np.arange(0,N,0.1),surfaceF(np.array(np.arange(0,N,0.1))))

plt.show()

