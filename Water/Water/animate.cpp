//
//  animate.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "animate.h"

#define SOLID 0
#define FLUID 1
#define EMPTY 2
int n_particles;
particle_t *particles;
float u[N+1][N+1][N+1],  v[N+1][N+1][N+1],  w[N+1][N+1][N+1],
      u2[N+1][N+1][N+1], v2[N+1][N+1][N+1], w2[N+1][N+1][N+1];
char types[N][N][N];
int d[N+1][N+1][N+1];
int fluid_cells_idx[N][N][N];
float dt = 0.2;
float gy = -0.4f;

void config1(void){
  n_particles = 0;
  for (int i=1; i<N/2; i++){
    for (int j=1; j<N/2; j++){
      for (int k=1; k<N/2; k++){
        types[i][j][k]=FLUID;
        n_particles+=8;
      }
    }
  }
}

void config2(void){
  n_particles = 0;
  /*for (int i=1; i<N-1; i++){
    for (int j=1; j<N/10; j++){
      for (int k=1; k<N-1; k++){
        types[i][j][k]=FLUID;
        n_particles+=8;
      }
    }
  }*/
  int cx=N/2;
  int cy=3*N/4;
  int cz=N/2;
  int r2=(N/8)*(N/8);
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N-1; j++){
      for (int k=1; k<N-1; k++){
        if ((i-cx)*(i-cx)+(j-cy)*(j-cy)+(k-cz)*(k-cz)<r2){
          types[i][j][k]=FLUID;
          n_particles+=8;
        }
      }
    }
  }
}

extern "C" void init_animation(void){
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N+1; k++){
        u[i][j][k]=0.0f;
        v[i][j][k]=0.0f;
        w[i][j][k]=0.0f;
      }
    }
  }
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        types[i][j][k]=EMPTY;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      types[0][i][j]=SOLID;
      types[i][0][j]=SOLID;
      types[i][j][0]=SOLID;
      types[N-1][i][j]=SOLID;
      types[i][N-1][j]=SOLID;
      types[i][j][N-1]=SOLID;
    }
  }
  config1();
  particles = (particle_t*)malloc(n_particles*sizeof(particle_t));
  int c = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        if (types[i][j][k]==FLUID){
          for (char d=0; d<8; d++){
            float dx=(d&1)*0.5+0.25;
            float dy=((d>>1)&1)*0.5+0.25;
            float dz=((d>>2)&1)*0.5+0.25;
            particles[c].x= i+dx+(rand2()*0.5-0.25);
            particles[c].y= j+dy+(rand2()*0.5-0.25);
            particles[c].z= k+dz+(rand2()*0.5-0.25);
            c++;
          }
        }
      }
    }
  }
}

float trilinear_interp(float x, float y, float z, float x0, float y0, float z0,
                       int i, int j, int k, float f[N+1][N+1][N+1]){
  float C000 = f[i][j][k];
  float C001 = f[i][j][k+1];
  float C010 = f[i][j+1][k];
  float C011 = f[i][j+1][k+1];
  float C100 = f[i+1][j][k];
  float C101 = f[i+1][j][k+1];
  float C110 = f[i+1][j+1][k];
  float C111 = f[i+1][j+1][k+1];
  float C00 = C000 * (x0+1-x) + C100 * (x-x0);
  float C01 = C001 * (x0+1-x) + C101 * (x-x0);
  float C10 = C010 * (x0+1-x) + C110 * (x-x0);
  float C11 = C011 * (x0+1-x) + C111 * (x-x0);
  float C0 = C00 * (y0+1-y) + C10 * (y-y0);
  float C1 = C01 * (y0+1-y) + C11 * (y-y0);
  return C0 * (z0+1-z) + C1 * (z-z0);
}
float trilinear_interp_u(float x, float y, float z){
  x=fminf(fmaxf(x,0),  N-0.000001);
  y=fminf(fmaxf(y,0.5),N-0.500001);
  z=fminf(fmaxf(z,0.5),N-0.500001);
  float x0 = (int)x;
  float y0 = (int)(y-0.5)+0.5;
  float z0 = (int)(z-0.5)+0.5;
  int i = (int)x;
  int j = (int)(y-0.5);
  int k = (int)(z-0.5);
  return trilinear_interp(x, y, z, x0, y0, z0, i, j, k, u);
}
float trilinear_interp_v(float x, float y, float z){
  x=fminf(fmaxf(x,0.5),N-0.500001);
  y=fminf(fmaxf(y,0),  N-0.000001);
  z=fminf(fmaxf(z,0.5),N-0.500001);
  float x0 = (int)(x-0.5)+0.5;
  float y0 = (int)y;
  float z0 = (int)(z-0.5)+0.5;
  int i = (int)(x-0.5);
  int j = (int)y;
  int k = (int)(z-0.5);
  return trilinear_interp(x, y, z, x0, y0, z0, i, j, k, v);
}
float trilinear_interp_w(float x, float y, float z){
  x=fminf(fmaxf(x,0.5),N-0.500001);
  y=fminf(fmaxf(y,0.5),N-0.500001);
  z=fminf(fmaxf(z,0),  N-0.000001);
  float x0 = (int)(x-0.5)+0.5;
  float y0 = (int)(y-0.5)+0.5;
  float z0 = (int)z;
  int i = (int)(x-0.5);
  int j = (int)(y-0.5);
  int k = (int)z;
  return trilinear_interp(x, y, z, x0, y0, z0, i, j, k, w);
}
// back is 1 if forward, -1 if backward
void RK2(float curr_x, float curr_y, float curr_z, float *ret_x, float *ret_y, float *ret_z, int back){
  float k1_x = trilinear_interp_u(curr_x, curr_y, curr_z);
  float k1_y = trilinear_interp_v(curr_x, curr_y, curr_z);
  float k1_z = trilinear_interp_w(curr_x, curr_y, curr_z);
  float temp_x = curr_x + 0.5*dt*k1_x*back;
  float temp_y = curr_y + 0.5*dt*k1_y*back;
  float temp_z = curr_z + 0.5*dt*k1_z*back;
  float k2_x = trilinear_interp_u(temp_x, temp_y, temp_z);
  float k2_y = trilinear_interp_v(temp_x, temp_y, temp_z);
  float k2_z = trilinear_interp_w(temp_x, temp_y, temp_z);
  *ret_x = curr_x + dt * k2_x*back;
  *ret_y = curr_y + dt * k2_y*back;
  *ret_z = curr_z + dt * k2_z*back;
}

void advect(void){
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        float curr_x = i;
        float curr_y = j+0.5;
        float curr_z = k+0.5;
        float prev_x, prev_y, prev_z;
        RK2(curr_x, curr_y, curr_z, &prev_x, &prev_y, &prev_z, -1);
        u2[i][j][k]=trilinear_interp_u(prev_x, prev_y, prev_z);
      }
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N; k++){
        float curr_x = i+0.5;
        float curr_y = j;
        float curr_z = k+0.5;
        float prev_x, prev_y, prev_z;
        RK2(curr_x, curr_y, curr_z, &prev_x, &prev_y, &prev_z, -1);
        v2[i][j][k]=trilinear_interp_v(prev_x, prev_y, prev_z);
      }
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N+1; k++){
        float curr_x = i+0.5;
        float curr_y = j+0.5;
        float curr_z = k;
        float prev_x, prev_y, prev_z;
        RK2(curr_x, curr_y, curr_z, &prev_x, &prev_y, &prev_z, -1);
        w2[i][j][k]=trilinear_interp_w(prev_x, prev_y, prev_z);
      }
    }
  }
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N+1; k++){
        u[i][j][k]=u2[i][j][k];
        v[i][j][k]=v2[i][j][k];
        w[i][j][k]=w2[i][j][k];
      }
    }
  }
}

void extrapolate(void){
  std::vector<int> W;
  int diff[] = {1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1};
  //U
  for (int i=0; i<N+1; i++)
    for (int j=0; j<N+1; j++)
      for (int k=0; k<N+1; k++)
        d[i][j][k]=INT_MAX;
  for (int i=0; i<N-1; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        if (types[i][j][k]==FLUID || types[i+1][j][k]==FLUID)
          d[i+1][j][k]=0;
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        if (d[i][j][k]==0) continue;
        for (int c=0; c<6; c++){
          int i2=i+diff[c*3];
          int j2=j+diff[c*3+1];
          int k2=k+diff[c*3+2];
          if (i2>=0 && i2<N+1 && j2>=0 && j2<N && k2>=0 && k2<N && d[i2][j2][k2]==0){
            d[i][j][k]=1;
            W.push_back(i);
            W.push_back(j);
            W.push_back(k);
            break;
          }
        }
      }
    }
  }
  int t=0;
  while (t<W.size()){
    int i=W[t];
    int j=W[t+1];
    int k=W[t+2];
    float sum = 0;
    int count = 0;
    for (int c=0; c<6; c++){
      int i2=i+diff[c*3];
      int j2=j+diff[c*3+1];
      int k2=k+diff[c*3+2];
      if (i2>=0 && i2<N+1 && j2>=0 && j2<N && k2>=0 && k2<N){
        if (d[i2][j2][k2]<d[i][j][k]){
          sum += u[i2][j2][k2];
          count++;
        } else if (d[i2][j2][k2]==INT_MAX){
          d[i2][j2][k2]=d[i][j][k]+1;
          W.push_back(i2);
          W.push_back(j2);
          W.push_back(k2);
        }
      }
    }
    u[i][j][k]=sum/count;
    t+=3;
  }
  W.clear();
  //V
  for (int i=0; i<N+1; i++)
    for (int j=0; j<N+1; j++)
      for (int k=0; k<N+1; k++)
        d[i][j][k]=INT_MAX;
  for (int i=0; i<N; i++)
    for (int j=0; j<N-1; j++)
      for (int k=0; k<N; k++)
      if (types[i][j][k]==FLUID || types[i][j+1][k]==FLUID)
        d[i][j+1][k]=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N; k++){
        if (d[i][j][k]==0) continue;
        for (int c=0; c<6; c++){
          int i2=i+diff[c*3];
          int j2=j+diff[c*3+1];
          int k2=k+diff[c*3+2];
          if (i2>=0 && i2<N && j2>=0 && j2<N+1 && k2>=0 && k2<N && d[i2][j2][k2]==0){
            d[i][j][k]=1;
            W.push_back(i);
            W.push_back(j);
            W.push_back(k);
            break;
          }
        }
      }
    }
  }
  t=0;
  while (t<W.size()){
    int i=W[t];
    int j=W[t+1];
    int k=W[t+2];
    float sum = 0;
    int count = 0;
    for (int c=0; c<6; c++){
      int i2=i+diff[c*3];
      int j2=j+diff[c*3+1];
      int k2=k+diff[c*3+2];
      if (i2>=0 && i2<N && j2>=0 && j2<N+1 && k2>=0 && k2<N){
        if (d[i2][j2][k2]<d[i][j][k]){
          sum += v[i2][j2][k2];
          count++;
        } else if (d[i2][j2][k2]==INT_MAX){
          d[i2][j2][k2]=d[i][j][k]+1;
          W.push_back(i2);
          W.push_back(j2);
          W.push_back(k2);
        }
      }
    }
    v[i][j][k]=sum/count;
    t+=3;
  }
  W.clear();
  //W
  for (int i=0; i<N+1; i++)
    for (int j=0; j<N+1; j++)
      for (int k=0; k<N+1; k++)
        d[i][j][k]=INT_MAX;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N-1; k++)
      if (types[i][j][k]==FLUID || types[i][j][k+1]==FLUID)
        d[i][j][k+1]=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N+1; k++){
        if (d[i][j][k]==0) continue;
        for (int c=0; c<6; c++){
          int i2=i+diff[c*3];
          int j2=j+diff[c*3+1];
          int k2=k+diff[c*3+2];
          if (i2>=0 && i2<N && j2>=0 && j2<N && k2>=0 && k2<N+1 && d[i2][j2][k2]==0){
            d[i][j][k]=1;
            W.push_back(i);
            W.push_back(j);
            W.push_back(k);
            break;
          }
        }
      }
    }
  }
  t=0;
  while (t<W.size()){
    int i=W[t];
    int j=W[t+1];
    int k=W[t+2];
    float sum = 0;
    int count = 0;
    for (int c=0; c<6; c++){
      int i2=i+diff[c*3];
      int j2=j+diff[c*3+1];
      int k2=k+diff[c*3+2];
      if (i2>=0 && i2<N && j2>=0 && j2<N && k2>=0 && k2<N+1){
        if (d[i2][j2][k2]<d[i][j][k]){
          sum += w[i2][j2][k2];
          count++;
        } else if (d[i2][j2][k2]==INT_MAX){
          d[i2][j2][k2]=d[i][j][k]+1;
          W.push_back(i2);
          W.push_back(j2);
          W.push_back(k2);
        }
      }
    }
    w[i][j][k]=sum/count;
    t+=3;
  }
}

void project(void){
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N+1; k++){
        u2[i][j][k]=0;
        v2[i][j][k]=0;
        w2[i][j][k]=0;
      }
    }
  }
        
  int n_fluid = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        if (types[i][j][k]==FLUID){
          n_fluid++;
        }
      }
    }
  }
  int *fluid_cells = (int *)malloc(n_fluid * 3 * sizeof(int));
  int idx=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        if (types[i][j][k]==FLUID){
          fluid_cells[idx*3]=i;
          fluid_cells[idx*3+1]=j;
          fluid_cells[idx*3+2]=k;
          fluid_cells_idx[i][j][k]=idx;
          idx++;
        }
      }
    }
  }
  mallocCG(n_fluid, 7);
  initCG();
  for (idx=0; idx<n_fluid; idx++){
    int i = fluid_cells[idx*3];
    int j = fluid_cells[idx*3+1];
    int k = fluid_cells[idx*3+2];
    if (types[i-1][j][k]==FLUID)
      write_A(idx, fluid_cells_idx[i-1][j][k], -1);
    if (types[i+1][j][k]==FLUID)
      write_A(idx, fluid_cells_idx[i+1][j][k], -1);
    if (types[i][j-1][k]==FLUID)
      write_A(idx, fluid_cells_idx[i][j-1][k], -1);
    if (types[i][j+1][k]==FLUID)
      write_A(idx, fluid_cells_idx[i][j+1][k], -1);
    if (types[i][j][k-1]==FLUID)
      write_A(idx, fluid_cells_idx[i][j][k-1], -1);
    if (types[i][j][k+1]==FLUID)
      write_A(idx, fluid_cells_idx[i][j][k+1], -1);
    int non_solid = 0;
    float D=0;
    if (types[i-1][j][k]!=SOLID){
      D-=u[i][j][k];
      non_solid++;
    }
    if (types[i+1][j][k]!=SOLID){
      D+=u[i+1][j][k];
      non_solid++;
    }
    if (types[i][j-1][k]!=SOLID){
      D-=v[i][j][k];
      non_solid++;
    }
    if (types[i][j+1][k]!=SOLID){
      D+=v[i][j+1][k];
      non_solid++;
    }
    if (types[i][j][k-1]!=SOLID){
      D-=w[i][j][k];
      non_solid++;
    }
    if (types[i][j][k+1]!=SOLID){
      D+=w[i][j][k+1];
      non_solid++;
    }
    write_A(idx, idx, non_solid);
    write_b(idx, -(1.0/dt)*D);
  }
  float *x = cg();
  for (int i=0; i<N-1; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        if (types[i][j][k]!=FLUID && types[i+1][j][k]!=FLUID)
          continue;
        float rightP=0;
        float leftP=0;
        if (types[i+1][j][k]==FLUID)
          rightP = x[fluid_cells_idx[i+1][j][k]];
        if (types[i+1][j][k]==EMPTY)
          rightP = 0;
        if (types[i+1][j][k]==SOLID){
          u2[i+1][j][k]=0;
          continue;
        }
        if (types[i][j][k]==FLUID)
          leftP = x[fluid_cells_idx[i][j][k]];
        if (types[i][j][k]==EMPTY)
          leftP = 0;
        if (types[i][j][k]==SOLID){
          u2[i+1][j][k] = 0;
          continue;
        }
        u2[i+1][j][k] = u[i+1][j][k] - dt * (rightP - leftP);
      }
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N-1; j++){
      for (int k=0; k<N; k++){
        if (types[i][j][k]!=FLUID && types[i][j+1][k]!=FLUID)
          continue;
        float topP = 0;
        float bottomP = 0;
        if (types[i][j+1][k]==FLUID)
          topP = x[fluid_cells_idx[i][j+1][k]];
        if (types[i][j+1][k]==EMPTY)
          topP = 0;
        if (types[i][j+1][k]==SOLID){
          v2[i][j+1][k] = 0;
          continue;
        }
        if (types[i][j][k]==FLUID)
          bottomP = x[fluid_cells_idx[i][j][k]];
        if (types[i][j][k]==EMPTY)
          bottomP = 0;
        if (types[i][j][k]==SOLID){
          v2[i][j+1][k] = 0;
          continue;
        }
        v2[i][j+1][k] = v[i][j+1][k] - dt * (topP - bottomP);
      }
    }
  }
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N-1; k++){
        if (types[i][j][k]!=FLUID && types[i][j][k+1]!=FLUID)
          continue;
        float frontP = 0;
        float backP = 0;
        if (types[i][j][k+1]==FLUID)
          frontP = x[fluid_cells_idx[i][j][k+1]];
        if (types[i][j][k+1]==EMPTY)
          frontP = 0;
        if (types[i][j][k+1]==SOLID){
          w2[i][j][k+1] = 0;
          continue;
        }
        if (types[i][j][k]==FLUID)
          backP = x[fluid_cells_idx[i][j][k]];
        if (types[i][j][k]==EMPTY)
          backP = 0;
        if (types[i][j][k]==SOLID){
          w2[i][j][k+1] = 0;
          continue;
        }
        w2[i][j][k+1] = w[i][j][k+1] - dt * (frontP - backP);
      }
    }
  }
  free(fluid_cells);
  freeCG();
  for (int i=0; i<N+1; i++){
    for (int j=0; j<N+1; j++){
      for (int k=0; k<N+1; k++){
        u[i][j][k]=u2[i][j][k];
        v[i][j][k]=v2[i][j][k];
        w[i][j][k]=w2[i][j][k];
      }
    }
  }
}

void move_particles(void){
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
        if (types[i][j][k]==FLUID)
          types[i][j][k]=EMPTY;
  for (int idx=0; idx<n_particles; idx++){
    float newX, newY, newZ;
    RK2(particles[idx].x,particles[idx].y,particles[idx].z,&newX,&newY,&newZ,1);
    if (types[(int)newX][(int)newY][(int)newZ]!=SOLID){
      particles[idx].x=newX;
      particles[idx].y=newY;
      particles[idx].z=newZ;
      types[(int)newX][(int)newY][(int)newZ]=FLUID;
    }
  }
}

extern "C" void animate(void){
  double t1 = timestamp();
  advect();
  project();
  extrapolate();
  move_particles();
  for (int i=0; i<N; i++)
    for (int j=0; j<N-1; j++)
      for (int k=0; k<N; k++)
        if (types[i][j][k]==FLUID || types[i][j+1][k]==FLUID)
          v[i][j+1][k] += dt * gy;
  double t2 = timestamp();
  printf("Animate: %f\n",t2-t1);
}
