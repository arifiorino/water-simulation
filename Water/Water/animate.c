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
float **u, **v, **u2, **v2;
char **types;
float *particles;
float dt = 0.2;
float gy = -0.4f;

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
  n_particles = 0;
  for (int i=1; i<N/2; i++){
    for (int j=1; j<N-1; j++){
      types[i][j]=FLUID;
      n_particles+=4;
    }
  }
  particles = (float*)malloc(n_particles*2*sizeof(float));
  int c = 0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==FLUID){
        particles[c] = i+0.25; particles[c+1] = j+0.25; c+=2;
        particles[c] = i+0.25; particles[c+1] = j+0.75; c+=2;
        particles[c] = i+0.75; particles[c+1] = j+0.25; c+=2;
        particles[c] = i+0.75; particles[c+1] = j+0.75; c+=2;
      }
    }
  }
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
int fluid_cell_idx(int i, int j, int *fluid_cells, int n_fluid){
  for (int idx=0; idx<n_fluid; idx++){
    if (fluid_cells[idx*2]==i && fluid_cells[idx*2+1]==j){
      return idx;
    }
  }
  return -1;
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
  int idx=0;
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==FLUID){
        fluid_cells[idx*2]=i;
        fluid_cells[idx*2+1]=j;
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
      write_A(idx, fluid_cell_idx(i-1,j,fluid_cells,n_fluid), -1);
    if (types[i][j-1]==FLUID)
      write_A(idx, fluid_cell_idx(i,j-1,fluid_cells,n_fluid), -1);
    if (types[i+1][j]==FLUID)
      write_A(idx, fluid_cell_idx(i+1,j,fluid_cells,n_fluid), -1);
    if (types[i][j+1]==FLUID)
      write_A(idx, fluid_cell_idx(i,j+1,fluid_cells,n_fluid), -1);
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
        rightP = x[fluid_cell_idx(i+1,j,fluid_cells,n_fluid)];
      if (types[i+1][j]==EMPTY)
        rightP = 0;
      if (types[i+1][j]==SOLID){
        u2[i+1][j]=0;
        continue;
      }
      if (types[i][j]==FLUID)
        leftP = x[fluid_cell_idx(i,j,fluid_cells,n_fluid)];
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
        topP = x[fluid_cell_idx(i,j+1,fluid_cells,n_fluid)];
      if (types[i][j+1]==EMPTY)
        topP = 0;
      if (types[i][j+1]==SOLID){
        v2[i][j+1] = 0;
        continue;
      }
      if (types[i][j]==FLUID)
        bottomP = x[fluid_cell_idx(i,j,fluid_cells,n_fluid)];
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
  freeCG();
}

void move_particles(void){
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if (types[i][j]==FLUID)
        types[i][j]=EMPTY;
  for (int idx=0; idx<n_particles; idx++){
    RK2(particles[idx*2],particles[idx*2+1],particles+idx*2,particles+idx*2+1);
    float i=particles[idx*2];
    float j=particles[idx*2+1];
    types[(int)i][(int)j]=FLUID;
  }
}

void animate(void){
  advect();
  project();
  move_particles();
  for (int i=0; i<N; i++)
    for (int j=0; j<N-1; j++)
      if (types[i][j]==FLUID || types[i][j+1]==FLUID)
        v[i][j+1] += dt * gy;
}

char water(float x, float y){
  //float F = 0;
  for (int i=0; i<n_particles; i++){
    float s = fabs(particles[i*2]-x) + fabs(particles[i*2+1]-y);
    if (s < 1) return 1;
      //F += (1-s*s)*(1-s*s)*(1-s*s);
  }
  //if (F<2) return 1;
  return 0;
}

void marching_squares(void){
  float gridSize = 0.5f;
  n_triangles = 0;
  for (float x=0; x<N-gridSize; x+=gridSize){
    for (float y=0; y<N-gridSize; y+=gridSize){
      char TL = water(x,y+gridSize);
      char TR = water(x+gridSize,y+gridSize);
      char BR = water(x+gridSize,y);
      char BL = water(x,y);
      if (TL + TR + BR + BL == 1)
        n_triangles++;
      else if (TL + TR + BR + BL == 2)
        n_triangles+=2;
      else if (TL + TR + BR + BL == 3)
        n_triangles+=3;
      else if (TL + TR + BR + BL == 4)
        n_triangles+=2;
    }
  }
  free(triangles);
  triangles = malloc(n_triangles * 3 * 2 * 4);
  int i=0;
  for (float x=0; x<N-gridSize; x+=gridSize){
    for (float y=0; y<N-gridSize; y+=gridSize){
      bool TL = water(x,y+gridSize)==1;
      bool TR = water(x+gridSize,y+gridSize)==1;
      bool BR = water(x+gridSize,y)==1;
      bool BL = water(x,y)==1;
      float x2 = x+gridSize;
      float y2 = y+gridSize;
      float xh = x+gridSize/2;
      float yh = y+gridSize/2;
      if       ( TL && TR && BR && BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=y;
        triangles[i+6]=x; triangles[i+7]=y; triangles[i+8]=x2; triangles[i+9]=y2; triangles[i+10]=x2; triangles[i+11]=y; i+=12;
      }else if ( TL &&!TR &&!BR &&!BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=xh; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=yh; i+=6;
      }else if (!TL && TR &&!BR &&!BL){
        triangles[i]=x2; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=yh; triangles[i+4]=xh; triangles[i+5]=y2; i+=6;
      }else if (!TL &&!TR && BR &&!BL){
        triangles[i]=x2; triangles[i+1]=y; triangles[i+2]=xh; triangles[i+3]=y; triangles[i+4]=x2; triangles[i+5]=yh; i+=6;
      }else if (!TL &&!TR &&!BR && BL){
        triangles[i]=x; triangles[i+1]=y; triangles[i+2]=x; triangles[i+3]=yh; triangles[i+4]=xh; triangles[i+5]=y; i+=6;
      }else if ( TL && TR &&!BR &&!BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=yh;
        triangles[i+6]=x; triangles[i+7]=yh; triangles[i+8]=x2; triangles[i+9]=y2; triangles[i+10]=x2; triangles[i+11]=yh; i+=12;
      }else if ( TL &&!TR && BR &&!BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=xh; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=yh;
        triangles[i+6]=x2; triangles[i+7]=y; triangles[i+8]=xh; triangles[i+9]=y; triangles[i+10]=x2; triangles[i+11]=yh; i+=12;
      }else if ( TL &&!TR &&!BR && BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=xh; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=y;
        triangles[i+6]=x; triangles[i+7]=y; triangles[i+8]=xh; triangles[i+9]=y2; triangles[i+10]=xh; triangles[i+11]=y; i+=12;
      }else if (!TL && TR && BR &&!BL){
        triangles[i]=xh; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=y2; triangles[i+4]=x2; triangles[i+5]=y;
        triangles[i+6]=xh; triangles[i+7]=y2; triangles[i+8]=x2; triangles[i+9]=y; triangles[i+10]=xh; triangles[i+11]=y; i+=12;
      }else if (!TL && TR &&!BR && BL){
        triangles[i]=xh; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=y2; triangles[i+4]=x2; triangles[i+5]=yh;
        triangles[i+6]=x; triangles[i+7]=yh; triangles[i+8]=xh; triangles[i+9]=y; triangles[i+10]=x; triangles[i+11]=y; i+=12;
      }else if (!TL &&!TR && BR && BL){
        triangles[i]=x; triangles[i+1]=y; triangles[i+2]=x; triangles[i+3]=yh; triangles[i+4]=x2; triangles[i+5]=yh;
        triangles[i+6]=x; triangles[i+7]=y; triangles[i+8]=x2; triangles[i+9]=yh; triangles[i+10]=x2; triangles[i+11]=y; i+=12;
      }else if (!TL && TR && BR && BL){
        triangles[i]=x; triangles[i+1]=y; triangles[i+2]=x; triangles[i+3]=yh; triangles[i+4]=x2; triangles[i+5]=y;
        triangles[i+6]=x; triangles[i+7]=yh; triangles[i+8]=xh; triangles[i+9]=y2; triangles[i+10]=x2; triangles[i+11]=y;
        triangles[i+12]=xh; triangles[i+13]=y2; triangles[i+14]=x2; triangles[i+15]=y2; triangles[i+16]=x2; triangles[i+17]=y; i+=18;
      }else if ( TL &&!TR && BR && BL){
        triangles[i]=x; triangles[i+1]=y; triangles[i+2]=x; triangles[i+3]=y2; triangles[i+4]=xh; triangles[i+5]=y2;
        triangles[i+6]=x; triangles[i+7]=y; triangles[i+8]=xh; triangles[i+9]=y2; triangles[i+10]=x2; triangles[i+11]=yh;
        triangles[i+12]=x; triangles[i+13]=y; triangles[i+14]=x2; triangles[i+15]=yh; triangles[i+16]=x2; triangles[i+17]=y; i+=18;
      }else if ( TL && TR &&!BR && BL){
        triangles[i]=x; triangles[i+1]=y; triangles[i+2]=x; triangles[i+3]=y2; triangles[i+4]=xh; triangles[i+5]=y;
        triangles[i+6]=x; triangles[i+7]=y2; triangles[i+8]=x2; triangles[i+9]=yh; triangles[i+10]=xh; triangles[i+11]=y;
        triangles[i+12]=x; triangles[i+13]=y2; triangles[i+14]=x2; triangles[i+15]=y2; triangles[i+16]=x2; triangles[i+17]=yh; i+=18;
      }else if ( TL && TR && BR &&!BL){
        triangles[i]=x; triangles[i+1]=y2; triangles[i+2]=x2; triangles[i+3]=y2; triangles[i+4]=x; triangles[i+5]=yh;
        triangles[i+6]=x; triangles[i+7]=yh; triangles[i+8]=x2; triangles[i+9]=y2; triangles[i+10]=xh; triangles[i+11]=y;
        triangles[i+12]=xh; triangles[i+13]=y; triangles[i+14]=x2; triangles[i+15]=y2; triangles[i+16]=x2; triangles[i+17]=y; i+=18;
      }
    }
  }
}
