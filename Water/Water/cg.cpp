//
//  cg.c
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include "cg.h"
float **A;
int **A_i;
int max_row;
float epsilon = 0.000001f;
float *b, *x;
float *scratch1, *scratch2, *scratch3, *scratch4;
int n;

void mallocCG(int n2, float max_row2){
  n = n2;
  max_row = max_row2;
  A = (float **)malloc2D(max_row, n, sizeof(float));
  A_i = (int **)malloc2D(max_row, n, sizeof(int));
  b = (float *)malloc(n * sizeof(float));
  x = (float *)malloc(n * sizeof(float));
  scratch1 = (float*)malloc(n * sizeof(float));
  scratch2 = (float*)malloc(n * sizeof(float));
  scratch3 = (float*)malloc(n * sizeof(float));
  scratch4 = (float*)malloc(n * sizeof(float));
}
void freeCG(void){
  free2D((void**)A, max_row);
  free2D((void**)A_i, max_row);
  free(b);
  free(x);
  free(scratch1);
  free(scratch2);
  free(scratch3);
  free(scratch4);
}

void initCG(){
  for(int i=0; i<max_row; i++){
    for(int j=0; j<n; j++){
      A_i[i][j]=-1;
    }
  }
  for(int i=0; i<n; i++){
    x[i]=0.0f;
  }
}
void write_A(int i, int j, float a){
  for (int k=0; k<max_row; k++){
    if (A_i[k][i] == -1){
      A_i[k][i] = j;
      A[k][i] = a;
      return;
    }
  }
  printf("Ran out!\n");
}
void write_b(int i, float a){
  b[i] = a;
}
// r = Av
void AMul(float *v, float *r){
  for (int i=0; i<n; i++){
    r[i]=0.0f;
    for (int k=0; k<max_row; k++){
      if (A_i[k][i]==-1) continue;
      r[i]+=A[k][i] * v[A_i[k][i]];
    }
  }
}
// r = u * v
float dot(float *u, float *v){
  float r = 0.0f;
  for(int i=0; i<n; i++){
    r+=u[i]*v[i];
  }
  return r;
}
// r = cv + u
void vecMulAdd(float *r, float c, float *v, float *u){
  for(int i=0; i<n; i++){
    r[i]=c*v[i]+u[i];
  }
}
//3461
float *cg(void){
  float *r = scratch1;
  float *p = scratch2;
  float *Ax = scratch3;
  float *Ap = scratch4;
  AMul(x, Ax);
  vecMulAdd(r, -1.0f, Ax, b);
  memcpy(p, r, n * sizeof(float));
  float rsold = dot(r, r);
  int i;
  for (i=0; i<n; i++){
    AMul(p, Ap);
    float alpha;
    if (rsold == 0.0f)
      alpha = 0.0f;
    else
      alpha = rsold / dot(p, Ap);
    vecMulAdd(x, alpha, p, x);
    vecMulAdd(r, -1.0f * alpha, Ap, r);
    float rsnew = dot(r, r);
    if (sqrt(rsnew) < epsilon){
      break;
    }
    vecMulAdd(p, rsnew / rsold, p, r);
    rsold = rsnew;
  }
  //printf("Iterations: %d\n",i);
  return x;
}

