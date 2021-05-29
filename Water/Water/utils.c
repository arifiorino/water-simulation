//
//  utils.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "utils.h"

void **malloc2D(int w, int h, int s){
  void **A = (void **)malloc(sizeof(void*)*w);
  for(int i=0; i < w; i++) {
    A[i] = (void *)malloc(h*s);
  }
  return A;
}
