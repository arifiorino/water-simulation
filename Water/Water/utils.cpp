//
//  utils.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "utils.h"
int n_threads = 8;
std::thread threads[7];

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

double rand2(void){
    return (double)rand() / (double)RAND_MAX ;
}

double timestamp(void){
  struct timespec spec;
  clock_gettime(CLOCK_REALTIME, &spec);
  return spec.tv_sec + spec.tv_nsec / 1e9;
}

void cross_prod(float ax, float ay, float az, float bx, float by, float bz,
           float *sx, float *sy, float *sz){
  *sx = ay * bz - az * by;
  *sy = az * bx - ax * bz;
  *sz = ax * by - ay * bx;
}

void parallel_for(int n, std::function<void (int start, int end)> functor){
  
  int batch_size = n/n_threads;
  for(int i = 0; i < n_threads-1; ++i){
    int start=i*batch_size;
    threads[i] = std::thread(functor, start, start+batch_size);
  }
  functor((n_threads-1)*batch_size, n);
  for(int i = 0; i < n_threads-1; ++i)
    threads[i].join();
}
