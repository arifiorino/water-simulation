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
void free2D(void **p, int w){
  for(int i=0; i < w; i++) {
    free(p[i]);
  }
  free(p);
}
void swap2D(float ***a, float ***b){
  float **tmp = *a;
  *a = *b;
  *b = tmp;
}
double rand2(void){
    return (double)rand() / (double)RAND_MAX ;
}

double timestamp(void){
  struct timespec spec;
  clock_gettime(CLOCK_REALTIME, &spec);
  return spec.tv_sec + spec.tv_nsec / 1e9;
}

void print2D(float **matrix, int w, int h){
  printf("[[");
  for (int i = 0; i < w; ++i){
    printf(" [");
    for (int j = 0; j < h; ++j)
      printf("%f, ", matrix[i][j]);
    printf("],\n");
  }
  printf("]\n");
}

void print2DInt(int **matrix, int w, int h){
  printf("[[");
  for (int i = 0; i < w; ++i){
    printf(" [");
    for (int j = 0; j < h; ++j)
      printf("%d, ", matrix[i][j]);
    printf("],\n");
  }
  printf("]\n");
}
