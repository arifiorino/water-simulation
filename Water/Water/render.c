//
//  render.c
//  Water
//
//  Created by Ari Fiorino on 6/5/21.
//

#include "render.h"
#include "animate.h"
#include "utils.h"

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

char water(float x, float y){
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
          if (F>0.3) return 1;
          curr = curr->next;
        }
      }
    }
  }
  return 0;
}

void marching_squares(void){
  hash_particles();
  float gridSize = 0.1;
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

void render(void){
  hash_particles();
  marching_squares();
}
