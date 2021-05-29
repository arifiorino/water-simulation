//
//  animate.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "animate.h"

#define BOUNDARY 0
#define FULL 1
#define SURFACE 2
#define EMPTY 3
float **p, **u, **v, **u2, **v2;
char **types;
float dt = 0.2;
float dx = 1;
float dy = 1;
float atmP = 0;
float waterP = 0.2;
float gx = 0;
float gy = 0;

void **malloc2D(int w, int h, int s){
  void **A = (void **)malloc(sizeof(void*)*w);
  for(int i=0; i < w; i++) {
    A[i] = (void *)malloc(s*h);
  }
  return A;
}

float interp(float **f, int w, int ax, int ay, int bx, int by){
  return (f[ax][ay]+f[bx][by])/2.0;
}

float interpU(float x, float y){
  float right=x+0.5;
  float left=x-0.5;
  float top=y+0.5;
  float bottom=y-0.5;
  float horizontal = (int)(y-0.5)+1.0;
  float vertical = ((int)x)+0.5;
  float TL = (vertical - left)*(top - horizontal);
  float TR = (right - vertical)*(top - horizontal);
  float BL = (vertical - left)*(horizontal - bottom);
  float BR = (right - vertical)*(horizontal - bottom);
  int i = (int)vertical;
  int j = (int)(horizontal-0.1);
  return TL * u[i][j+1] + TR * u[i+1][j+1] + BL * u[i][j] + BR * u[i+1][j];
}

float interpV(float x, float y){
  float right=x+0.5;
  float left=x-0.5;
  float top=y+0.5;
  float bottom=y-0.5;
  float horizontal = ((int)y)+0.5;
  float vertical = ((int)x-0.5)+1.0;
  float TL = (vertical - left)*(top - horizontal);
  float TR = (right - vertical)*(top - horizontal);
  float BL = (vertical - left)*(horizontal - bottom);
  float BR = (right - vertical)*(horizontal - bottom);
  int i = (int)(vertical-0.1);
  int j = (int)horizontal;
  return TL * v[i][j+1] + TR * v[i+1][j+1] + BL * v[i][j] + BR * v[i+1][j];
}

void setBoundarySurface(void){
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      types[i][j] = EMPTY;
    }
  }
  for (int i=0; i<NParticles; i++){
    int x = (int)(particles[i*2]);
    int y = (int)(particles[i*2+1]);
    types[x][y] = FULL;
  }
  for (int i=0; i<N; i++){
    types[0][i]=BOUNDARY;
    types[i][0]=BOUNDARY;
    types[N-1][i]=BOUNDARY;
    types[i][N-1]=BOUNDARY;
  }
  for (int i=1; i<N-1; i++){
    for (int j=1; j<N-1; j++){
      bool top = types[i][j+1] == EMPTY;
      bool right = types[i+1][j] == EMPTY;
      bool bottom = types[i][j-1] == EMPTY;
      bool left = types[i-1][j] == EMPTY;
      if (types[i][j]==FULL && (left || right || top || bottom)){
        types[i][j]=SURFACE;
      }
    }
  }
  //Set border vel - normal
  for (int i=1; i<N-1; i++){
    u[0][i] = -u[2][i];
    u[N][i] = -u[N-2][i];
    v[i][0] = -v[i][2];
    v[i][N] = -v[i][N-2];
  }
  //Set border vel - wall
  for (int i=1; i<N-1; i++){
    u[1][i] = 0;
    u[N-1][i] = 0;
    v[i][1] = 0;
    v[i][N-1] = 0;
  }
  //Set border vel - tangential
  for (int i=1; i<N; i++){
    u[i][0] = u[i][1];
    u[i][N-1] = u[i][N-2];
    v[0][i] = v[1][i];
    v[N-1][i] = v[N-2][i];
  }
  u[0][0]=0;
  v[0][0]=0;
  u[N][0]=0;
  v[0][N]=0;
  u[0][N-1]=0;
  v[N-1][0]=0;
  u[N][N-1]=0;
  v[N-1][N]=0;

  //Set border pressure
  for (int i=1; i<N-1; i++){
    p[0][i]=p[1][i];
    p[i][0]=p[i][1];
    p[N-1][i]=p[N-2][i];
    p[i][N-1]=p[i][N-2];
  }
  //Set empty pressure, vel
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      if (types[i][j]==EMPTY){
        u[i+1][j]=0;
        u[i][j]=0;
        v[i][j+1]=0;
        v[i][j]=0;
        p[i][j]=atmP;
      }
    }
  }
}

void setSurface(void){
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      bool L = types[i-1][j]==EMPTY;
      bool B = types[i][j-1]==EMPTY;
      bool T = types[i][j+1]==EMPTY;
      bool R = types[i+1][j]==EMPTY;
      if      ( L&&!B&&!T&&!R){ v2[i-1][j] = -v[i][j]-v[i][j+1]+u[i+1][j];}
      else if (!L&& B&&!T&&!R){ u2[i][j-1] = -u[i][j]+v[i][j+1]+u[i+1][j];}
      else if (!L&&!B&& T&&!R){ v2[i][j+1] = u[i][j]+v[i][j]-u[i+1][j];}
      else if (!L&&!B&&!T&& R){ u2[i+1][j] = u[i][j]+v[i][j]-v[i][j+1];}

      else if ( L&& B&&!T&&!R){ u2[i][j] = u[i+1][j]; v2[i][j]=v[i][j+1];}
      else if ( L&&!B&& T&&!R){ u2[i][j] = u[i+1][j]; v2[i][j+1]=v[i][j];}
      else if ( L&&!B&&!T&& R){}
      else if (!L&& B&& T&&!R){}
      else if (!L&& B&&!T&& R){ v2[i][j]=v[i][j+1]; u2[i+1][j]=u[i][j];}
      else if (!L&&!B&& T&& R){ v2[i][j+1]=v[i][j]; u2[i+1][j]=u[i][j];}

      else if (!L&&!B&&!T&& R){ u2[i][j]=u[i+1][j];}
      else if (!L&&!B&& T&&!R){ v2[i][j]=v[i][j+1];}
      else if (!L&& B&&!T&&!R){ v2[i][j+1]=v[i][j];}
      else if ( L&&!B&&!T&&!R){ u2[i+1][j]=u[i][j];}

      else if ( L&& B&& T&& R){ u2[i][j]=gx;v2[i][j]=gy;v2[i][j+1]=gy;u2[i+1][j]=gx;}
      p[i][j]=atmP;
    }
  }
  //Switch u-u2, v-v2
  float **tmp = u2;
  u2 = u;
  u = tmp;
  tmp = v2;
  v2 = v;
  v = tmp;
}

void pressureSolve(void){
  
}


void initAnimation(void){
  p = (float**)malloc2D(N, N, sizeof(float));
  u = (float**)malloc2D(N+1, N, sizeof(float));
  u2 = (float**)malloc2D(N+1, N, sizeof(float));
  v = (float**)malloc2D(N, N+1, sizeof(float));
  v2 = (float**)malloc2D(N, N+1, sizeof(float));
  types = (char**)malloc2D(N, N, sizeof(char));
  NParticles = 72;
  particles = (float*)malloc(NParticles*2*sizeof(float));
}

void animate(void){
  
}