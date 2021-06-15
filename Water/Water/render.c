//
//  render.c
//  Water
//
//  Created by Ari Fiorino on 6/5/21.
//

#include "render.h"
#include "animate.h"
#include "utils.h"

#define split 5
bool level_set[N*split+1][N*split+1][N*split+1];
int indices_size, vertices_size;

typedef struct vertex_ll{
  float x; float y; float z;
  struct vertex_ll *next;
} vertex_ll_t;
bool vertices_hash[N*split+1][N*split+1][N*split+1];

void init_render(void){
  n_vertices = 0;
  vertices_size = 0;
  n_indices = 0;
  indices_size = 0;
  indices = NULL;
  vertices = NULL;
  normals = NULL;
}

char tetrahedra_to_triangles[16][12] =
 {{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  { 0, 1, 0, 3, 0, 2,-1,-1,-1,-1,-1,-1},
  { 0, 1, 1, 2, 1, 3,-1,-1,-1,-1,-1,-1},
  { 0, 2, 1, 3, 0, 3, 0, 2, 1, 2, 1, 3},
  { 0, 2, 2, 3, 1, 2,-1,-1,-1,-1,-1,-1},
  { 0, 1, 0, 3, 2, 3, 0, 1, 2, 3, 1, 2},
  { 0, 1, 0, 2, 2, 3, 0, 1, 2, 3, 1, 3},
  { 0, 3, 2, 3, 1, 3,-1,-1,-1,-1,-1,-1},
  { 0, 3, 1, 3, 2, 3,-1,-1,-1,-1,-1,-1},
  { 0, 1, 2, 3, 0, 2, 0, 1, 1, 3, 2, 3},
  { 0, 1, 2, 3, 0, 3, 0, 1, 1, 2, 2, 3},
  { 0, 2, 1, 2, 2, 3,-1,-1,-1,-1,-1,-1},
  { 0, 2, 0, 3, 1, 3, 0, 2, 1, 3, 1, 2},
  { 0, 1, 1, 3, 1, 2,-1,-1,-1,-1,-1,-1},
  { 0, 1, 0, 2, 0, 3,-1,-1,-1,-1,-1,-1},
  {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};

char cube_to_tetrahedra[2][5][12] =
 {{{0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1},
   {0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0},
   {1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0},
   {1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0},
   {1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0}},
  {{0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1},
   {1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0},
   {1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0},
   {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0},
   {1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0}}};

void hash_particles(void){
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      for (int k=0; k<N; k++)
      particles_hash[i][j][k]=NULL;
  for (int idx=0; idx<n_particles; idx++){
    int i=particles[idx].x;
    int j=particles[idx].y;
    int k=particles[idx].z;
    particles[idx].next=NULL;
    if (particles_hash[i][j][k] == NULL){
      particles_hash[i][j][k]=&(particles[idx]);
      continue;
    }
    particle_t *curr = particles_hash[i][j][k];
    while (curr->next!=NULL) curr=curr->next;
    curr->next = particles+idx;
  }
}

bool metaballs(float x, float y, float z){
  float F = 0;
  for (int di=-1; di<=1; di++){
    for (int dj=-1; dj<=1; dj++){
      for (int dk=-1; dk<=1; dk++){
        int i=x; int j=y; int k=z;
        if (i+di>=0 && i+di<N && j+dj>=0 && j+dj<N && k+dk>=0 && k+dk<N){
          particle_t *curr = particles_hash[i+di][j+dj][k+dk];
          while (curr != NULL){
            float s2 = (curr->x-x)*(curr->x-x) + (curr->y-y)*(curr->y-y) + (curr->z-z)*(curr->z-z);
            if (s2 < 1) F += (1-s2)*(1-s2)*(1-s2);
            if (F > 0.8) return true;
            curr = curr->next;
          }
        }
      }
    }
  }
  return false;
}

void add_vertex(float x, float y, float z){
  
}

/*
void add_triangle(float tri[9]){
  n_triangles++;
  if (n_triangles*3*3 >= triangles_size){
    triangles_size *= 2;
    triangles = realloc(triangles, triangles_size * sizeof(float));
    printf("%p\n",triangles);
  }
  for (int i=0; i<9; i++)
    triangles[(n_triangles-1)*3*3+i]=tri[i] * gridSize;
}

void tetrahedra(int tetra[12]){
  bool a = points[tetra[0]][tetra[1]][tetra[2]];
  bool b = points[tetra[3]][tetra[4]][tetra[5]];
  bool c = points[tetra[6]][tetra[7]][tetra[8]];
  bool d = points[tetra[9]][tetra[10]][tetra[11]];
  int idx = d + (c<<1) + (b<<2) + (a<<3);
  for (int idx2=0; idx2<=6 && lookup[idx][idx2] != -1; idx2+=6){
    float tri[9];
    for (int t = 0; t < 3; t++){
      char p0 = lookup[idx][idx2+t*2];
      char p1 = lookup[idx][idx2+t*2+1];
      float x = 0;
      float y = 0;
      float z = 0;
      if (p0==0 || p1==0){ x+=tetra[0]; y+=tetra[1]; z+=tetra[2]; }
      if (p0==1 || p1==1){ x+=tetra[3]; y+=tetra[4]; z+=tetra[5]; }
      if (p0==2 || p1==2){ x+=tetra[6]; y+=tetra[7]; z+=tetra[8]; }
      if (p0==3 || p1==3){ x+=tetra[9]; y+=tetra[10]; z+=tetra[11]; }
      x/=2; y/=2; z/=2;
      tri[t*3]=x; tri[t*3+1]=y; tri[t*3+2]=z;
    }
    add_triangle(tri);
  }
}

void marching_tetrahedra(void){
  n_triangles = 0;
  int count = 0;
  for (int i=0; i<N*gridSplit+1; i++){
    for (int j=0; j<N*gridSplit+1; j++){
      for (int k=0; k<N*gridSplit+1; k++){
      points[i][j][k]=water(i*gridSize,j*gridSize,k*gridSize);
      count++;
      }
    }
  }
  for (int i=0; i<N*gridSplit; i++){
    for (int j=0; j<N*gridSplit; j++){
      for (int k=0; k<N*gridSplit; k++){
        for (int tetra_i = 0; tetra_i<5; tetra_i++){
          int tetra[12];
          for (int q=0; q<12; q++)
            tetra[q]=tetrahedras[(i+j+k)%2][tetra_i][q];
          tetra[0]+=i; tetra[1]+=j; tetra[2]+=k;
          tetra[3]+=i; tetra[4]+=j; tetra[5]+=k;
          tetra[6]+=i; tetra[7]+=j; tetra[8]+=k;
          tetra[9]+=i; tetra[10]+=j; tetra[11]+=k;
          tetrahedra(tetra);
        }
      }
    }
  }
}
*/

void render(void){
  hash_particles();
  double t1 = timestamp();
  //marching_tetrahedra();
  double t2 = timestamp();
  printf("March:   %f\n",t2-t1);
}
