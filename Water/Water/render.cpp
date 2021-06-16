//
//  render.c
//  Water
//
//  Created by Ari Fiorino on 6/5/21.
//

#include "render.h"

#define split 5
bool level_set[N*split+1][N*split+1][N*split+1];
int indices_size, vertices_size;
std::unordered_map<int, int> vertex_to_index;
std::vector <float> vertices_v, normals_v;
std::vector<int> indices_v;

int n_indices, n_vertices;
int *indices;
float *vertices, *normals;

int max_vertex_i = (N * split + 1) * 2;

char tetrahedron_to_triangles[16][12] =
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

char cube_to_tetrahedron[2][5][4][3] =
 {{{{0, 1, 1}, {1, 0, 1}, {0, 0, 0}, {0, 0, 1}},
   {{0, 0, 0}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}},
   {{1, 0, 1}, {0, 1, 1}, {0, 0, 0}, {1, 1, 0}},
   {{1, 0, 0}, {1, 0, 1}, {0, 0, 0}, {1, 1, 0}},
   {{1, 1, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}}},
  {{{0, 0, 1}, {1, 1, 1}, {1, 0, 0}, {1, 0, 1}},
   {{1, 1, 1}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}},
   {{1, 1, 1}, {0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
   {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0}},
   {{1, 1, 1}, {0, 1, 1}, {0, 0, 1}, {0, 1, 0}}}};

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
  return (x-8)*(x-8) + (y-8)*(y-8) + (z-8)*(z-8) < 16;
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

int hash_vertex(int vertex[3]){
  return vertex[0] + (vertex[1] * max_vertex_i) +
         (vertex[2] * max_vertex_i * max_vertex_i);
}

int add_vertex(int vertex[3]){
  int h = hash_vertex(vertex);
  try{
    return vertex_to_index.at(h);
  }catch (const std::out_of_range& oor){
    vertices_v.push_back((float)vertex[0]/split/2.0);
    vertices_v.push_back((float)vertex[1]/split/2.0);
    vertices_v.push_back((float)vertex[2]/split/2.0);
    normals_v.push_back(0);
    normals_v.push_back(0);
    normals_v.push_back(0);
    int index = (int)(vertices_v.size()/3)-1;
    vertex_to_index.insert({h, index});
    return index;
  }
}

void add_triangle(int triangle[3][3]){
  int index1 = add_vertex(triangle[0]);
  int index2 = add_vertex(triangle[1]);
  int index3 = add_vertex(triangle[2]);
  float ax = vertices_v[index1*3];
  float ay = vertices_v[index1*3+1];
  float az = vertices_v[index1*3+2];
  float bx = vertices_v[index2*3];
  float by = vertices_v[index2*3+1];
  float bz = vertices_v[index2*3+2];
  float cx = vertices_v[index3*3];
  float cy = vertices_v[index3*3+1];
  float cz = vertices_v[index3*3+2];
  float nx, ny, nz;
  cross_prod(bx-ax, by-ay, bz-az, cx-bx, cy-by, cz-bz, &nx, &ny, &nz);
  indices_v.push_back(index1);
  indices_v.push_back(index2);
  indices_v.push_back(index3);
  normals_v[index1*3]+=nx; normals_v[index1*3+1]+=ny; normals_v[index1*3+2]+=nz;
  normals_v[index2*3]+=nx; normals_v[index2*3+1]+=ny; normals_v[index2*3+2]+=nz;
  normals_v[index3*3]+=nx; normals_v[index3*3+1]+=ny; normals_v[index3*3+2]+=nz;
}

void add_tetrahedron(int tetrahedron[4][3]){
  bool b1 = level_set[tetrahedron[0][0]][tetrahedron[0][1]][tetrahedron[0][2]];
  bool b2 = level_set[tetrahedron[1][0]][tetrahedron[1][1]][tetrahedron[1][2]];
  bool b3 = level_set[tetrahedron[2][0]][tetrahedron[2][1]][tetrahedron[2][2]];
  bool b4 = level_set[tetrahedron[3][0]][tetrahedron[3][1]][tetrahedron[3][2]];
  int tetrahedron_idx = b1 + (b2<<1) + (b3<<2) + (b4<<3);
  for (int triangle_i=0; triangle_i<12 && tetrahedron_to_triangles[tetrahedron_idx][triangle_i]!=-1;
       triangle_i+=6){
    int triangle[3][3];
    for (int t = 0; t < 3; t++){
      char p0 = tetrahedron_to_triangles[tetrahedron_idx][triangle_i + t*2];
      char p1 = tetrahedron_to_triangles[tetrahedron_idx][triangle_i + t*2+1];
      int x = 0; int y = 0; int z = 0;
      if (p0==0 || p1==0){ x+=tetrahedron[0][0]; y+=tetrahedron[0][1]; z+=tetrahedron[0][2]; }
      if (p0==1 || p1==1){ x+=tetrahedron[1][0]; y+=tetrahedron[1][1]; z+=tetrahedron[1][2]; }
      if (p0==2 || p1==2){ x+=tetrahedron[2][0]; y+=tetrahedron[2][1]; z+=tetrahedron[2][2]; }
      if (p0==3 || p1==3){ x+=tetrahedron[3][0]; y+=tetrahedron[3][1]; z+=tetrahedron[3][2]; }
      triangle[t][0]=x; triangle[t][1]=y; triangle[t][2]=z;
    }
    add_triangle(triangle);
  }
}

void add_cube(int cube[3]){
  for (int tetrahedron_i = 0; tetrahedron_i<5; tetrahedron_i++){
    int tetrahedron[4][3];
    for (int i=0; i<4; i++)
      for (int j=0; j<3; j++)
        tetrahedron[i][j]=cube_to_tetrahedron[(cube[0]+cube[1]+cube[2])%2][tetrahedron_i][i][j];
    tetrahedron[0][0]+=cube[0]; tetrahedron[0][1]+=cube[1]; tetrahedron[0][2]+=cube[2];
    tetrahedron[1][0]+=cube[0]; tetrahedron[1][1]+=cube[1]; tetrahedron[1][2]+=cube[2];
    tetrahedron[2][0]+=cube[0]; tetrahedron[2][1]+=cube[1]; tetrahedron[2][2]+=cube[2];
    tetrahedron[3][0]+=cube[0]; tetrahedron[3][1]+=cube[1]; tetrahedron[3][2]+=cube[2];
    add_tetrahedron(tetrahedron);
  }
}


void marching_tetrahedra(void){
  for (int i=0; i<N*split+1; i++){
    for (int j=0; j<N*split+1; j++){
      for (int k=0; k<N*split+1; k++){
      level_set[i][j][k]=metaballs((float)i/split,(float)j/split,(float)k/split);
      }
    }
  }
  for (int i=0; i<N*split; i++){
    for (int j=0; j<N*split; j++){
      for (int k=0; k<N*split; k++){
        int cube[] = {i,j,k};
        add_cube(cube);
      }
    }
  }
  
}

extern "C" void render(void){
  vertex_to_index.clear();
  vertices_v.clear();
  normals_v.clear();
  indices_v.clear();
  
  double t1 = timestamp();
  hash_particles();
  marching_tetrahedra();
  double t2 = timestamp();
  n_vertices = (int)(vertices_v.size()/3);
  n_indices = (int)indices_v.size();
  vertices = vertices_v.data();
  normals = normals_v.data();
  indices = indices_v.data();
  printf("Render:    %f\n",t2-t1);
}
