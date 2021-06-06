//
//  animate.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "animate.h"
#include "utils.h"
#include "cg.h"

#define SOLID 0
#define FLUID 1
#define EMPTY 2
int N = 64;
int n_particles;
particle_t *particles;
particle_t ***particles_hash;
float **u, **v, **u2, **v2;
char **types;
float dt = 0.2;
float gy = -0.4f;

void config1(void){
  n_particles = 0;
  for (int i=1; i<N/2; i++){
    for (int j=1; j<N-1; j++){
      types[i][j]=FLUID;
      n_particles+=4;
    }
  }
}

void config2(void){
  n_particles = 0;
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N/8; j++){
      types[i][j]=FLUID;
      n_particles+=4;
    }
  }
  int cx=N/2;
  int cy=3*N/4;
  int r2=(N/8)*(N/8);
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N-1; j++){
      if ((i-cx)*(i-cx)+(j-cy)*(j-cy)<r2){
        types[i][j]=FLUID;
        n_particles+=4;
      }
    }
  }
}

void config3(void){
  n_particles = 0;
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N/4; j++){
      types[i][j]=FLUID;
      n_particles+=4;
    }
  }
  int cx1=N/4;
  int cx2=N/2;
  int cx3=3*N/4;
  int cy=3*N/4;
  int r2=(N/8)*(N/8);
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N-1; j++){
      if ((i-cx1)*(i-cx1)+(j-cy)*(j-cy)<r2 ||
          (i-cx2)*(i-cx2)+(j-cy)*(j-cy)<r2 ||
          (i-cx3)*(i-cx3)+(j-cy)*(j-cy)<r2){
        types[i][j]=FLUID;
        n_particles+=4;
      }
    }
  }
}

void initAnimation(void){
  u = (float**)malloc2D(N+1, N, sizeof(float));
  u2 = (float**)malloc2D(N+1, N, sizeof(float));
  v = (float**)malloc2D(N, N+1, sizeof(float));
  v2 = (float**)malloc2D(N, N+1, sizeof(float));
  types = (char**)malloc2D(N, N, sizeof(char));
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N; j++){
      u[i][j]=0.0f;
      v[j][i]=0.0f;
    }
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      types[i][j]=EMPTY;
  for (int i=0; i<N; i++){
    types[0][i]=SOLID;
    types[i][0]=SOLID;
    types[N-1][i]=SOLID;
    types[i][N-1]=SOLID;
  }
  config1();
  particles = (particle_t*)malloc(n_particles*sizeof(particle_t));
  int c = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==FLUID){
        particles[c].x= i+0.25+(rand2()*0.5-0.25); particles[c].y= j+0.25+(rand2()*0.5-0.25); c++;
        particles[c].x= i+0.25+(rand2()*0.5-0.25); particles[c].y= j+0.75+(rand2()*0.5-0.25); c++;
        particles[c].x= i+0.75+(rand2()*0.5-0.25); particles[c].y= j+0.25+(rand2()*0.5-0.25); c++;
        particles[c].x= i+0.75+(rand2()*0.5-0.25); particles[c].y= j+0.75+(rand2()*0.5-0.25); c++;
      }
    }
  }
  particles_hash = (particle_t***)malloc2D(N, N, sizeof(particle_t*));
}
float bilinear_interp_u(float x, float y){
  x=fmaxf(x,0);
  x=fminf(x,N-0.0001);
  y=fmaxf(y,0.5);
  y=fminf(y,N-0.5001);
  float x1 = (float)(int)x;
  float x2 = x1+1;
  float y1 = (int)(y-0.5)+0.5;
  float y2 = y1+1;
  int i = (int)x;
  int j = (int)(y-0.5);
  float Q11 = u[i][j];
  float Q12 = u[i][j+1];
  float Q21 = u[i+1][j];
  float Q22 = u[i+1][j+1];
  return Q11*(x2-x)*(y2-y) + Q21*(x-x1)*(y2-y) +
         Q12*(x2-x)*(y-y1) + Q22*(x-x1)*(y-y1);
}
float bilinear_interp_v(float x, float y){
  x=fmaxf(x,0.5);
  x=fminf(x,N-0.5001);
  y=fmaxf(y,0);
  y=fminf(y,N-0.0001);
  float x1 = (int)(x-0.5)+0.5;
  float x2 = x1+1;
  float y1 = (int)y;
  float y2 = y1+1;
  int i = (int)(x-0.5);
  int j = (int)y;
  float Q11 = v[i][j];
  float Q12 = v[i][j+1];
  float Q21 = v[i+1][j];
  float Q22 = v[i+1][j+1];
  return Q11*(x2-x)*(y2-y) + Q21*(x-x1)*(y2-y) +
         Q12*(x2-x)*(y-y1) + Q22*(x-x1)*(y-y1);
}
void RK2(float curr_x, float curr_y, float *ret_x, float *ret_y){
  float k1_x = bilinear_interp_u(curr_x, curr_y);
  float k1_y = bilinear_interp_v(curr_x, curr_y);
  float temp_x = curr_x + 0.5*dt*k1_x;
  float temp_y = curr_y + 0.5*dt*k1_y;
  float k2_x = bilinear_interp_u(temp_x, temp_y);
  float k2_y = bilinear_interp_v(temp_x, temp_y);
  *ret_x = curr_x + dt * k2_x;
  *ret_y = curr_y + dt * k2_y;
}
void backwards_RK2(float curr_x, float curr_y, float *ret_x, float *ret_y){
  float k1_x = bilinear_interp_u(curr_x, curr_y);
  float k1_y = bilinear_interp_v(curr_x, curr_y);
  float temp_x = curr_x - 0.5*dt*k1_x;
  float temp_y = curr_y - 0.5*dt*k1_y;
  float k2_x = bilinear_interp_u(temp_x, temp_y);
  float k2_y = bilinear_interp_v(temp_x, temp_y);
  *ret_x = curr_x - dt * k2_x;
  *ret_y = curr_y - dt * k2_y;
}
void advect(void){
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N; j++){
      float curr_x = i;
      float curr_y = j+0.5;
      float prev_x;
      float prev_y;
      backwards_RK2(curr_x, curr_y, &prev_x, &prev_y);
      u2[i][j]=bilinear_interp_u(prev_x, prev_y);
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N+1; j++){
      float curr_x = i+0.5;
      float curr_y = j;
      float prev_x;
      float prev_y;
      backwards_RK2(curr_x, curr_y, &prev_x, &prev_y);
      v2[i][j]=bilinear_interp_v(prev_x, prev_y);
    }
  }
  swap2D(&u, &u2);
  swap2D(&v, &v2);
}

void extrapolate(void){
  int **d = (int **)malloc2D(N+1,N,sizeof(int));
  for (int i=0; i<N+1; i++)
    for (int j=0; j<N; j++)
      d[i][j]=INT_MAX;
  int *W = malloc(N*N*2*sizeof(int));
  int W_len=0;
  for (int i=0; i<N-1; i++)
    for (int j=0; j<N; j++)
      if (types[i][j]==FLUID || types[i+1][j]==FLUID)
        d[i+1][j]=0;
  int diff[] = {-1,0,0,-1,1,0,0,1};
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N; j++){
      if (d[i][j]==0) continue;
      for (int k=0; k<4; k++){
        int i2=i+diff[k*2];
        int j2=j+diff[k*2+1];
        if (i2>=0 && i2<N+1 && j2>=0 && j2<N && d[i2][j2]==0){
          d[i][j]=1;
          W[W_len]=i;
          W[W_len+1]=j;
          W_len+=2;
          break;
        }
      }
    }
  }
  int t=0;
  while (t<W_len){
    int i=W[t];
    int j=W[t+1];
    float new = 0;
    int count = 0;
    for (int k=0; k<4; k++){
      int i2=i+diff[k*2];
      int j2=j+diff[k*2+1];
      if (i2>=0 && i2<N+1 && j2>=0 && j2<N){
        if (d[i2][j2]<d[i][j]){
          new += u[i2][j2];
          count++;
        } else if (d[i2][j2]==INT_MAX){
          d[i2][j2]=d[i][j]+1;
          W[W_len]=i2;
          W[W_len+1]=j2;
          W_len+=2;
        }
      }
    }
    u[i][j]=new/count;
    t+=2;
  }
  free2D((void **)d, N+1);
  d = (int **)malloc2D(N,N+1,sizeof(int));
  for (int i=0; i<N; i++)
    for (int j=0; j<N+1; j++)
      d[i][j]=INT_MAX;
  for (int i=0; i<N; i++)
    for (int j=0; j<N-1; j++)
      if (types[i][j]==FLUID || types[i][j+1]==FLUID)
        d[i][j+1]=0;
  W_len=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N+1; j++){
      if (d[i][j]==0) continue;
      for (int k=0; k<4; k++){
        int i2=i+diff[k*2];
        int j2=j+diff[k*2+1];
        if (i2>=0 && i2<N && j2>=0 && j2<N+1 && d[i2][j2]==0){
          d[i][j]=1;
          W[W_len]=i;
          W[W_len+1]=j;
          W_len+=2;
          break;
        }
      }
    }
  }
  t=0;
  while (t<W_len){
    int i=W[t];
    int j=W[t+1];
    float new = 0;
    int count = 0;
    for (int k=0; k<4; k++){
      int i2=i+diff[k*2];
      int j2=j+diff[k*2+1];
      if (i2>=0 && i2<N && j2>=0 && j2<N+1){
        if (d[i2][j2]<d[i][j]){
          new += v[i2][j2];
          count++;
        } else if (d[i2][j2]==INT_MAX){
          d[i2][j2]=d[i][j]+1;
          W[W_len]=i2;
          W[W_len+1]=j2;
          W_len+=2;
        }
      }
    }
    v[i][j]=new/count;
    t+=2;
  }
  free2D((void **)d, N);
  free(W);
  
}

void project(void){
  int n_fluid = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==FLUID){
        n_fluid++;
      }
    }
  }
  int *fluid_cells = malloc(n_fluid * 2 * sizeof(int));
  int **fluid_cells_idx = (int **)malloc2D(N, N, sizeof(int));
  int idx=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==FLUID){
        fluid_cells[idx*2]=i;
        fluid_cells[idx*2+1]=j;
        fluid_cells_idx[i][j]=idx;
        idx++;
      }
    }
  }
  mallocCG(n_fluid, 5);
  initCG();
  for (idx=0; idx<n_fluid; idx++){
    int i = fluid_cells[idx*2];
    int j = fluid_cells[idx*2+1];
    if (types[i-1][j]==FLUID)
      write_A(idx, fluid_cells_idx[i-1][j], -1);
    if (types[i][j-1]==FLUID)
      write_A(idx, fluid_cells_idx[i][j-1], -1);
    if (types[i+1][j]==FLUID)
      write_A(idx, fluid_cells_idx[i+1][j], -1);
    if (types[i][j+1]==FLUID)
      write_A(idx, fluid_cells_idx[i][j+1], -1);
    int non_solid = 0;
    float D=0;
    if (types[i-1][j]!=SOLID){
      D-=u[i][j];
      non_solid+=1;
    }
    if (types[i][j-1]!=SOLID){
      D-=v[i][j];
      non_solid+=1;
    }
    if (types[i+1][j]!=SOLID){
      D+=u[i+1][j];
      non_solid+=1;
    }
    if (types[i][j+1]!=SOLID){
      D+=v[i][j+1];
      non_solid+=1;
    }
    write_A(idx, idx, non_solid);
    write_b(idx, -(1.0/dt)*D);
  }
  float *x = cg();
  for (int i=0; i<N-1; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]!=FLUID && types[i+1][j]!=FLUID)
        continue;
      float rightP=0;
      float leftP=0;
      if (types[i+1][j]==FLUID)
        rightP = x[fluid_cells_idx[i+1][j]];
      if (types[i+1][j]==EMPTY)
        rightP = 0;
      if (types[i+1][j]==SOLID){
        u2[i+1][j]=0;
        continue;
      }
      if (types[i][j]==FLUID)
        leftP = x[fluid_cells_idx[i][j]];
      if (types[i][j]==EMPTY)
        leftP = 0;
      if (types[i][j]==SOLID){
        u2[i+1][j] = 0;
        continue;
      }
      u2[i+1][j] = u[i+1][j] - dt * (rightP - leftP);
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N-1; j++){
      if (types[i][j]!=FLUID && types[i][j+1]!=FLUID)
        continue;
      float topP = 0;
      float bottomP = 0;
      if (types[i][j+1]==FLUID)
        topP = x[fluid_cells_idx[i][j+1]];
      if (types[i][j+1]==EMPTY)
        topP = 0;
      if (types[i][j+1]==SOLID){
        v2[i][j+1] = 0;
        continue;
      }
      if (types[i][j]==FLUID)
        bottomP = x[fluid_cells_idx[i][j]];
      if (types[i][j]==EMPTY)
        bottomP = 0;
      if (types[i][j]==SOLID){
        v2[i][j+1] = 0;
        continue;
      }
      v2[i][j+1] = v[i][j+1] - dt * (topP - bottomP);
    }
  }
  swap2D(&u, &u2);
  swap2D(&v, &v2);
  free(fluid_cells);
  free2D((void **)fluid_cells_idx, N);
  freeCG();
}

void move_particles(void){
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if (types[i][j]==FLUID)
        types[i][j]=EMPTY;
  for (int idx=0; idx<n_particles; idx++){
    float newX, newY;
    RK2(particles[idx].x,particles[idx].y,&newX,&newY);
    if (types[(int)newX][(int)newY]!=SOLID){
      particles[idx].x=newX;
      particles[idx].y=newY;
      types[(int)newX][(int)newY]=FLUID;
    }
  }
}

void animate(void){
  advect();
  double t1 = timestamp();
  project();
  double t2 = timestamp();
  extrapolate();
  move_particles();
  printf("Project: %f\n",t2-t1);
  for (int i=0; i<N; i++)
    for (int j=0; j<N-1; j++)
      if (types[i][j]==FLUID || types[i][j+1]==FLUID)
        v[i][j+1] += dt * gy;
}
