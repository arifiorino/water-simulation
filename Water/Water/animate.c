//
//  animate.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "animate.h"

void initAnimation(void){
  vertices = malloc(6*4);
  vertices[0] = -0.8f;
  vertices[1] = -0.8f;
  vertices[2] = 0.0f;
  vertices[3] = 0.8f;
  vertices[4] = 0.8f;
  vertices[5] = -0.8f;
}

void animate(void){
  vertices[0]+=0.001;
}
