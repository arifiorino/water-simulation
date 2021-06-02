import matplotlib.pyplot as plt
import numpy as np
import time

N = 128

#Calculate actual phi
phi = np.zeros((N, N))
surfaceF = lambda x: x/4+(np.sin(2*np.pi*x/(N/2)))*N/8+N/4
surface = [(x, surfaceF(x)) for x in np.arange(0,N,0.3)]
for i in range(N):
  for j in range(N):
    p=(i+0.5, j+0.5)
    phi[i][j]=np.sqrt(min([(p[0]-surface[k][0])**2 + (p[1]-surface[k][1])**2 for k in range(len(surface))]))
    if surfaceF(i+0.5) > j+0.5:
      phi[i][j]*=-1

start=time.time()
phi3 = np.zeros((N, N))
for i in range(N):
  for j in range(N):
    p=(i+0.5, j+0.5)
    v=np.sqrt(min([(p[0]-surface[k][0])**2 + (p[1]-surface[k][1])**2 for k in range(len(surface))]))
    phi3[i][j]=v
end=time.time()
print('simple',end-start)

start=time.time()
phi2 = np.zeros((N, N))
findZero = lambda y1, y2: -y1 / (y2-y1)
surface=[]
closest_point = [[None for i in range(N)] for j in range(N)]
known = [[False for i in range(N)] for j in range(N)]
for i in range(N):
  for j in range(N):
    distances=[]
    candidates=[]
    if i<N-1 and (phi[i][j]<0) != (phi[i+1][j]<0):
      d=findZero(phi[i][j],phi[i+1][j])
      distances.append(d)
      candidates.append((i+0.5+d,j+0.5))
    if i>0 and (phi[i][j]<0) != (phi[i-1][j]<0):
      d=findZero(phi[i][j],phi[i-1][j])
      distances.append(d)
      candidates.append((i+0.5-d,j+0.5))
    if j<N-1 and (phi[i][j]<0) != (phi[i][j+1]<0):
      d=findZero(phi[i][j],phi[i][j+1])
      distances.append(d)
      candidates.append((i+0.5,j+0.5+d))
    if j>0 and (phi[i][j]<0) != (phi[i][j-1]<0):
      d=findZero(phi[i][j],phi[i][j-1])
      distances.append(d)
      candidates.append((i+0.5,j+0.5-d))
    if len(candidates)>0:
      known[i][j]=True
      closest_point[i][j]=candidates[np.argmin(distances)]
      surface.append(closest_point[i][j])
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
end=time.time()
print('time',end-start)

diffs=[]
for i in range(N):
  for j in range(N):
    diffs.append(abs(phi2[i][j]-abs(phi[i][j])))
print('diff',sum(diffs)/len(diffs))
plt.pcolormesh(phi2.T)
plt.colorbar()
plt.plot(np.arange(0,N,0.1),surfaceF(np.array(np.arange(0,N,0.1))))
plt.show()

