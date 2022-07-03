//
//  render.c
//  Water
//
//  Created by Ari Fiorino on 6/5/21.
//

#include "render.h"
#include "animate.h"
#include "utils.h"

float *triangles;
int n_triangles;

void hash_particles(void){
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      particles_hash[i][j]=NULL;
  for (int idx=0; idx<n_particles; idx++){
    int i=particles[idx].x;
    int j=particles[idx].y;
    particles[idx].next=NULL;
    if (particles_hash[i][j] == NULL){
      particles_hash[i][j]=&(particles[idx]);
      continue;
    }
    particle_t *curr = particles_hash[i][j];
    while (curr->next!=NULL) curr=curr->next;
    curr->next = particles+idx;
  }
}

bool water(float x, float y){
  float F = 0;
  for (int di=-1; di<=1; di++){
    for (int dj=-1; dj<=1; dj++){
      int i=x;
      int j=y;
      if (i+di>=0 && i+di<N && j+dj>=0 && j+dj<N){
        particle_t *curr = particles_hash[i+di][j+dj];
        while (curr != NULL){
          float s = (curr->x-x)*(curr->x-x) + (curr->y-y)*(curr->y-y);
          if (s < 1) F += (1-s*s)*(1-s*s)*(1-s*s);
          if (F>0.8) return true;
          curr = curr->next;
        }
      }
    }
  }
  return false;
}

void marching_squares(void){
  hash_particles();
  int gridSplit = 5;
  float gridSize = 1.0f/gridSplit;
  n_triangles = 0;
  bool **points = (bool **)malloc2D(N*gridSplit+1, N*gridSplit+1, sizeof(bool));
  for (int i=0; i<N*gridSplit+1; i++){
    for (int j=0; j<N*gridSplit+1; j++){
      points[i][j]=water(i*gridSize,j*gridSize);
    }
  }
  for (int i=0; i<N*gridSplit; i++){
    for (int j=0; j<N*gridSplit; j++){
      bool TL = points[i][j+1];
      bool TR = points[i+1][j+1];
      bool BR = points[i+1][j];
      bool BL = points[i][j];
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
  triangles = (float*)malloc(n_triangles * 3 * 2 * 4);
  int c=0;
  for (int i=0; i<N*gridSplit; i++){
    for (int j=0; j<N*gridSplit; j++){
      bool TL = points[i][j+1];
      bool TR = points[i+1][j+1];
      bool BR = points[i+1][j];
      bool BL = points[i][j];
      float x = i*gridSize;
      float y = j*gridSize;
      float xh = x+gridSize/2;
      float yh = y+gridSize/2;
      float x2 = x+gridSize;
      float y2 = y+gridSize;
      x=(x/N/gridSplit/gridSize*2)-1;
      y=(y/N/gridSplit/gridSize*2)-1;
      xh=(xh/N/gridSplit/gridSize*2)-1;
      yh=(yh/N/gridSplit/gridSize*2)-1;
      x2=(x2/N/gridSplit/gridSize*2)-1;
      y2=(y2/N/gridSplit/gridSize*2)-1;
      
      if       ( TL && TR && BR && BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=y;
        triangles[c+6]=x; triangles[c+7]=y; triangles[c+8]=x2; triangles[c+9]=y2; triangles[c+10]=x2; triangles[c+11]=y; c+=12;
      }else if ( TL &&!TR &&!BR &&!BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=xh; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=yh; c+=6;
      }else if (!TL && TR &&!BR &&!BL){
        triangles[c]=x2; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=yh; triangles[c+4]=xh; triangles[c+5]=y2; c+=6;
      }else if (!TL &&!TR && BR &&!BL){
        triangles[c]=x2; triangles[c+1]=y; triangles[c+2]=xh; triangles[c+3]=y; triangles[c+4]=x2; triangles[c+5]=yh; c+=6;
      }else if (!TL &&!TR &&!BR && BL){
        triangles[c]=x; triangles[c+1]=y; triangles[c+2]=x; triangles[c+3]=yh; triangles[c+4]=xh; triangles[c+5]=y; c+=6;
      }else if ( TL && TR &&!BR &&!BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=yh;
        triangles[c+6]=x; triangles[c+7]=yh; triangles[c+8]=x2; triangles[c+9]=y2; triangles[c+10]=x2; triangles[c+11]=yh; c+=12;
      }else if ( TL &&!TR && BR &&!BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=xh; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=yh;
        triangles[c+6]=x2; triangles[c+7]=y; triangles[c+8]=xh; triangles[c+9]=y; triangles[c+10]=x2; triangles[c+11]=yh; c+=12;
      }else if ( TL &&!TR &&!BR && BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=xh; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=y;
        triangles[c+6]=x; triangles[c+7]=y; triangles[c+8]=xh; triangles[c+9]=y2; triangles[c+10]=xh; triangles[c+11]=y; c+=12;
      }else if (!TL && TR && BR &&!BL){
        triangles[c]=xh; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=y2; triangles[c+4]=x2; triangles[c+5]=y;
        triangles[c+6]=xh; triangles[c+7]=y2; triangles[c+8]=x2; triangles[c+9]=y; triangles[c+10]=xh; triangles[c+11]=y; c+=12;
      }else if (!TL && TR &&!BR && BL){
        triangles[c]=xh; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=y2; triangles[c+4]=x2; triangles[c+5]=yh;
        triangles[c+6]=x; triangles[c+7]=yh; triangles[c+8]=xh; triangles[c+9]=y; triangles[c+10]=x; triangles[c+11]=y; c+=12;
      }else if (!TL &&!TR && BR && BL){
        triangles[c]=x; triangles[c+1]=y; triangles[c+2]=x; triangles[c+3]=yh; triangles[c+4]=x2; triangles[c+5]=yh;
        triangles[c+6]=x; triangles[c+7]=y; triangles[c+8]=x2; triangles[c+9]=yh; triangles[c+10]=x2; triangles[c+11]=y; c+=12;
      }else if (!TL && TR && BR && BL){
        triangles[c]=x; triangles[c+1]=y; triangles[c+2]=x; triangles[c+3]=yh; triangles[c+4]=x2; triangles[c+5]=y;
        triangles[c+6]=x; triangles[c+7]=yh; triangles[c+8]=xh; triangles[c+9]=y2; triangles[c+10]=x2; triangles[c+11]=y;
        triangles[c+12]=xh; triangles[c+13]=y2; triangles[c+14]=x2; triangles[c+15]=y2; triangles[c+16]=x2; triangles[c+17]=y; c+=18;
      }else if ( TL &&!TR && BR && BL){
        triangles[c]=x; triangles[c+1]=y; triangles[c+2]=x; triangles[c+3]=y2; triangles[c+4]=xh; triangles[c+5]=y2;
        triangles[c+6]=x; triangles[c+7]=y; triangles[c+8]=xh; triangles[c+9]=y2; triangles[c+10]=x2; triangles[c+11]=yh;
        triangles[c+12]=x; triangles[c+13]=y; triangles[c+14]=x2; triangles[c+15]=yh; triangles[c+16]=x2; triangles[c+17]=y; c+=18;
      }else if ( TL && TR &&!BR && BL){
        triangles[c]=x; triangles[c+1]=y; triangles[c+2]=x; triangles[c+3]=y2; triangles[c+4]=xh; triangles[c+5]=y;
        triangles[c+6]=x; triangles[c+7]=y2; triangles[c+8]=x2; triangles[c+9]=yh; triangles[c+10]=xh; triangles[c+11]=y;
        triangles[c+12]=x; triangles[c+13]=y2; triangles[c+14]=x2; triangles[c+15]=y2; triangles[c+16]=x2; triangles[c+17]=yh; c+=18;
      }else if ( TL && TR && BR &&!BL){
        triangles[c]=x; triangles[c+1]=y2; triangles[c+2]=x2; triangles[c+3]=y2; triangles[c+4]=x; triangles[c+5]=yh;
        triangles[c+6]=x; triangles[c+7]=yh; triangles[c+8]=x2; triangles[c+9]=y2; triangles[c+10]=xh; triangles[c+11]=y;
        triangles[c+12]=xh; triangles[c+13]=y; triangles[c+14]=x2; triangles[c+15]=y2; triangles[c+16]=x2; triangles[c+17]=y; c+=18;
      }
    }
  }
  free2D((void **)points, N*gridSplit+1);
}

void render(void){
  hash_particles();
  double t1 = timestamp();
  marching_squares();
  double t2 = timestamp();
  printf("March:   %f\n",t2-t1);
}
