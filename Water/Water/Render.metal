//
//  Render.metal
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include <metal_stdlib>
using namespace metal;

#define N 32
#define split 8
#define scanN 512
#define scanlgN 9
#define sphereN 2103

inline void scan_2(threadgroup int* A, threadgroup int* B, int j, constant int *powLookup){
  threadgroup_barrier(mem_flags::mem_threadgroup);
  int x = 0;
  int y = 0;
  for (int i=0; i<scanlgN; i++){
    if (j >= powLookup[i]){
      x = A[j - powLookup[i]];
      y = B[j - powLookup[i]];
      threadgroup_barrier(mem_flags::mem_threadgroup);
      A[j] += x;
      B[j] += y;
      threadgroup_barrier(mem_flags::mem_threadgroup);
    }else{
      threadgroup_barrier(mem_flags::mem_threadgroup);
      threadgroup_barrier(mem_flags::mem_threadgroup);
    }
  }
}

inline int metaballs(float x1, float y1, float z1, float x2, float y2, float z2){
  float s2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
  if (s2 < 1)
    return (int)((1-s2)*(1-s2)*(1-s2)*100000);
  return 0;
}

kernel void calc_level_set(device float* particles [[buffer(0)]],
                           device atomic_int* level_set [[buffer(1)]],
                           constant int* sphere [[buffer(2)]],
                           uint3 pos [[thread_position_in_grid]]){
  int idx = pos.x;
  int i=particles[idx*3]*split;
  int j=particles[idx*3+1]*split;
  int k=particles[idx*3+2]*split;
  for (int sphereI=0; sphereI<sphereN*3; sphereI+=3){
    int di=sphere[sphereI];
    int dj=sphere[sphereI+1];
    int dk=sphere[sphereI+2];
    const int M = N * split + 1;
    if (i+di>=0 && k+dk>=0 && j+dj>=0 && i+di<M && j+dj<M && k+dk<M){
      int x =(i+di)*M*M+(j+dj)*M+(k+dk);
      int z =atomic_load_explicit(level_set+x, memory_order_relaxed);
      if (z>80000)
        continue;
      int m = metaballs(particles[idx*3], particles[idx*3+1], particles[idx*3+2],
                          (float)(i+di)/split, (float)(j+dj)/split, (float)(k+dk)/split);
      atomic_fetch_add_explicit(level_set+x, m, memory_order_relaxed);
    }
  }
}

inline int vertex_idx(int v[3]){
  return v[2]*(N*split*2)*(N*split*2) + v[1]*(N*split*2) + v[0];
}
inline int cube_idx_2(int c[3]){
  return c[2]*(N*split+1)*(N*split+1) + c[1]*(N*split+1) + c[0];
}
inline int cube_idx(int c[3]){
  return c[2]*(N*split)*(N*split) + c[1]*(N*split) + c[0];
}
inline int cube_to_tetrahedron_idx(int a, int b, int c, int d){
  return a*3*4*5+b*3*4+c*3+d;
}

inline void add_vertex(int v[3], int c[3], device int *vertex_to_index, thread int *n_vertices){
  if (c[0] == v[0]/2 && c[1] == v[1]/2 && c[2] == v[2]/2){
    int vi = vertex_idx(v);
    if (vertex_to_index[vi]==-1){
        vertex_to_index[vi]=*n_vertices;
        *n_vertices+=1;
    }
  }
}

inline void add_triangle(int triangle[3][3], int c[3], device int *vertex_to_index, thread int *n_vertices,
                         thread int *n_indices){
  add_vertex(triangle[0], c, vertex_to_index, n_vertices);
  add_vertex(triangle[1], c, vertex_to_index, n_vertices);
  add_vertex(triangle[2], c, vertex_to_index, n_vertices);
  *n_indices+=3;
}

inline void add_tetrahedron(int tetrahedron[4][3],
                            int c[3],
                            device int *level_set,
                            constant int *tetrahedron_to_triangles,
                            device int *vertex_to_index,
                            thread int *n_vertices,
                            thread int *n_indices,
                            device int *debug_arr){
  int a1 = cube_idx_2(tetrahedron[0]);
  int a2 = cube_idx_2(tetrahedron[1]);
  int a3 = cube_idx_2(tetrahedron[2]);
  int a4 = cube_idx_2(tetrahedron[3]);
  bool b1 = level_set[a1]>80000;
  bool b2 = level_set[a2]>80000;
  bool b3 = level_set[a3]>80000;
  bool b4 = level_set[a4]>80000;
  int tetrahedron_idx = ((int)b1) + ((int)b2<<1) + ((int)b3<<2) + ((int)b4<<3);
  
  for (int triangle_i=0; triangle_i<12; triangle_i+=6){
    if (tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i]!=-1){
      int triangle[3][3];
      for (int t = 0; t < 3; t++){
        char p0 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2];
        char p1 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2+1];
        int x = 0; int y = 0; int z = 0;
        if (p0==0 || p1==0){ x+=tetrahedron[0][0]; y+=tetrahedron[0][1]; z+=tetrahedron[0][2]; }
        if (p0==1 || p1==1){ x+=tetrahedron[1][0]; y+=tetrahedron[1][1]; z+=tetrahedron[1][2]; }
        if (p0==2 || p1==2){ x+=tetrahedron[2][0]; y+=tetrahedron[2][1]; z+=tetrahedron[2][2]; }
        if (p0==3 || p1==3){ x+=tetrahedron[3][0]; y+=tetrahedron[3][1]; z+=tetrahedron[3][2]; }
        triangle[t][0]=x; triangle[t][1]=y; triangle[t][2]=z;
      }
      add_triangle(triangle, c, vertex_to_index, n_vertices, n_indices);
    }
  }
}

//Threadgroup: 32kb
kernel void add_cube(threadgroup int* tg_scan [[threadgroup(0)]],
                     threadgroup int* tg_scan_2 [[threadgroup(1)]],
                     uint3 thread_position_in_grid [[thread_position_in_grid]],
                     uint3 thread_position_in_threadgroup [[thread_position_in_threadgroup]],
                     device int* level_set [[buffer(0)]],
                     device int* vertex_to_index[[buffer(1)]],
                     constant int* tetrahedron_to_triangles[[buffer(2)]],
                     constant int* cube_to_tetrahedron[[buffer(3)]],
                     device int* debug_arr [[buffer(4)]],
                     constant int* pow_lookup[[buffer(5)]],
                     device int* indices_scan[[buffer(6)]]){
  int n_vertices=0;
  int n_indices=0;
  int c[3];
  c[0] = thread_position_in_grid.x;
  c[1] = thread_position_in_grid.y;
  c[2] = thread_position_in_grid.z;
  for (int di=0; di<2; di++){
    for (int dj=0; dj<2; dj++){
      for (int dk=0; dk<2; dk++){
        int v[3];
        v[0]=c[0]*2+di;
        v[1]=c[1]*2+dj;
        v[2]=c[2]*2+dk;
        vertex_to_index[vertex_idx(v)] = -1;
      }
    }
  }
  
  for (int tetrahedron_i = 0; tetrahedron_i<5; tetrahedron_i++){
    int tetrahedron[4][3];
    for (int i=0; i<4; i++){
      for (int j=0; j<3; j++){
        tetrahedron[i][j]=cube_to_tetrahedron[cube_to_tetrahedron_idx((c[0]+c[1]+c[2])%2,tetrahedron_i,i,j)];
      }
    }
    tetrahedron[0][0]+=c[0]; tetrahedron[0][1]+=c[1]; tetrahedron[0][2]+=c[2];
    tetrahedron[1][0]+=c[0]; tetrahedron[1][1]+=c[1]; tetrahedron[1][2]+=c[2];
    tetrahedron[2][0]+=c[0]; tetrahedron[2][1]+=c[1]; tetrahedron[2][2]+=c[2];
    tetrahedron[3][0]+=c[0]; tetrahedron[3][1]+=c[1]; tetrahedron[3][2]+=c[2];
    add_tetrahedron(tetrahedron, c, level_set, tetrahedron_to_triangles, vertex_to_index, &n_vertices, &n_indices, debug_arr);
  }
  int g[3];
  g[0] = c[0]%8;
  g[1] = c[1]%8;
  g[2] = c[2]%8;
  int idx = g[0] + g[1]*8 + g[2]*64;
  int v[3];
  v[0] = c[0]*2;
  v[1] = c[1]*2;
  v[2] = c[2]*2;
  tg_scan[idx]=n_vertices;
  tg_scan_2[idx]=n_indices;
  scan_2(tg_scan, tg_scan_2, idx, pow_lookup);
  vertex_to_index[vertex_idx(v)]=tg_scan[idx];
  indices_scan[cube_idx(c)] = tg_scan_2[idx];
}


kernel void globalScan(device int* vertex_to_index[[buffer(0)]],
                       device int* indices_scan[[buffer(1)]]){
  int n_verts = 0;
  for (int k=14; k<N*split*2; k+=16){
    for (int j=14; j<N*split*2; j+=16){
      for (int i=14; i<N*split*2; i+=16){
        int v[3];
        v[0]=i; v[1]=j; v[2]=k;
        int idx = vertex_idx(v);
        n_verts += vertex_to_index[idx];
        vertex_to_index[idx] = n_verts;
      }
    }
  }
  int n_idx = 0;
  for (int k=7; k<N*split; k+=8){
    for (int j=7; j<N*split; j+=8){
      for (int i=7; i<N*split; i+=8){
        int c[3];
        c[0]=i; c[1]=j; c[2]=k;
        int idx = cube_idx(c);
        n_idx += indices_scan[idx];
        indices_scan[idx] = n_idx;
      }
    }
  }
}

//PASS 2!!!!!!!!!


int get_vertex_idx(int v[3], device int* vertex_to_index, thread int *a1, thread int *a2, thread int *a3,
                   thread int *a4, thread int *a5, thread int *a6, thread int *a7){
  int c[3];
  c[0]=v[0]/2;
  c[1]=v[1]/2;
  c[2]=v[2]/2;
  int box[3];
  box[0]=c[0]/8;
  box[1]=c[1]/8;
  box[2]=c[2]/8;
  int prevBox = 0;
  if (box[0]!=0 || box[1]!=0 || box[2]!=0){
    int idx = box[2]*32*32 + box[1]*32 + box[0] - 1;
    int pos[3];
    pos[0] = 2*((idx % 32) * 8 + 7);
    pos[1] = 2*(((idx/32)%32) * 8 + 7);
    pos[2] = 2*((idx/32/32) * 8 + 7);
    int x = pos[2]*(N*split*2)*(N*split*2) + pos[1]*(N*split*2) + pos[0];
    *a1 = pos[0];
    *a2 = pos[1];
    *a3 = pos[2];
    *a4 = x;
    prevBox = vertex_to_index[x];
  }
  int grid[3];
  grid[0] = c[0]%8;
  grid[1] = c[1]%8;
  grid[2] = c[2]%8;
  int prevGrid = 0;
  if (grid[0]!=0 || grid[1]!=0 || grid[2]!=0){
    int idx = grid[2]*8*8 + grid[1]*8 + grid[0] - 1;
    int pos[3];
    pos[0] = 2*(box[0]*8 + (idx % 8));
    pos[1] = 2*(box[1]*8 + ((idx/8)%8));
    pos[2] = 2*(box[2]*8 + (idx/8/8));
    int x = pos[2]*(N*split*2)*(N*split*2) + pos[1]*(N*split*2) + pos[0];
    prevGrid = vertex_to_index[x];
  }
  int x = v[2]*(N*split*2)*(N*split*2) + v[1]*(N*split*2) + v[0];
  int vi = vertex_to_index[x];
  *a5 = prevBox;
  *a6 = prevGrid;
  *a7 = vi;
  return prevBox + prevGrid + vi;
}

int get_indices_idx(int c[3], device int* indices_scan){
  int box[3];
  box[0]=c[0]/8;
  box[1]=c[1]/8;
  box[2]=c[2]/8;
  int prevBox = 0;
  if (box[0]!=0 || box[1]!=0 || box[2]!=0){
    int idx = box[2]*32*32 + box[1]*32 + box[0] - 1;
    int pos[3];
    pos[0] = (idx % 32) * 8 + 7;
    pos[1] = ((idx/32)%32) * 8 + 7;
    pos[2] = (idx/32/32) * 8 + 7;
    int x = pos[2]*(N*split)*(N*split) + pos[1]*(N*split) + pos[0];
    prevBox = indices_scan[x];
  }
  int grid[3];
  grid[0] = c[0]%8;
  grid[1] = c[1]%8;
  grid[2] = c[2]%8;
  int prevGrid = 0;
  if (grid[0]!=0 || grid[1]!=0 || grid[2]!=0){
    int idx = grid[2]*8*8 + grid[1]*8 + grid[0] - 1;
    int pos[3];
    pos[0] = box[0]*8 + (idx % 8);
    pos[1] = box[1]*8 + ((idx/8)%8);
    pos[2] = box[2]*8 + (idx/8/8);
    int x = pos[2]*(N*split)*(N*split) + pos[1]*(N*split) + pos[0];
    prevGrid = indices_scan[x];
  }
  return prevBox + prevGrid;
}

inline void cross_prod(int ax, int ay, int az, int bx, int by, int bz,
                       thread int *sx, thread int *sy, thread int *sz){
  *sx = ay * bz - az * by;
  *sy = az * bx - ax * bz;
  *sz = ax * by - ay * bx;
}

inline int add_vertex_2(int v[3], int c[3], device int* vertex_to_index,
                        device int* vertices, device int *debug_arr){
  int a1, a2, a3, a4, a5, a6, a7;
  int idx = get_vertex_idx(v, vertex_to_index, &a1, &a2, &a3, &a4, &a5 ,&a6, &a7);
  //debug_arr[cube_idx(c)] = -1;
  if (c[0] == v[0]/2 && c[1] == v[1]/2 && c[2] == v[2]/2){ //It is your vertex
    if (idx < 300000){ //3305016
      debug_arr[0] = v[0];
      debug_arr[1] = v[1];
      debug_arr[2] = v[2];
      debug_arr[3] = c[0];
      debug_arr[4] = c[1];
      debug_arr[5] = c[2];
      debug_arr[6] = idx;
      debug_arr[7]=a1;
      debug_arr[8]=a2;
      debug_arr[9]=a3;
      debug_arr[10]=a4;
      debug_arr[11]=a5;
      debug_arr[12]=a6;
      debug_arr[13]=a7;
    }
    vertices[idx*3  ]=v[0];
    vertices[idx*3+1]=v[1];
    vertices[idx*3+2]=v[2];
  }
  return idx;
}

inline void add_triangle_2(int triangle[3][3], device int *vertex_to_index,
                           device int *vertices,
                           device atomic_int *normals,
                           device int *indices,
                           thread int *index_idx, int c[3],device int* debug_arr){
  int index1 = add_vertex_2(triangle[0], c, vertex_to_index, vertices, debug_arr);
  int index2 = add_vertex_2(triangle[1], c, vertex_to_index, vertices, debug_arr);
  int index3 = add_vertex_2(triangle[2], c, vertex_to_index, vertices, debug_arr);
  indices[*index_idx] = index1;
  indices[*index_idx + 1] = index2;
  indices[*index_idx + 2] = index3;
  *index_idx+=3;
  int ax = vertices[index1*3];
  int ay = vertices[index1*3+1];
  int az = vertices[index1*3+2];
  int bx = vertices[index2*3];
  int by = vertices[index2*3+1];
  int bz = vertices[index2*3+2];
  int cx = vertices[index3*3];
  int cy = vertices[index3*3+1];
  int cz = vertices[index3*3+2];
  int nx, ny, nz;
  cross_prod(bx-ax, by-ay, bz-az, cx-bx, cy-by, cz-bz, &nx, &ny, &nz);
  atomic_fetch_add_explicit(normals+index1*3  , nx, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index1*3+1, ny, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index1*3+2, nz, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index2*3  , nx, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index2*3+1, ny, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index2*3+2, nz, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index3*3  , nx, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index3*3+1, ny, memory_order_relaxed);
  atomic_fetch_add_explicit(normals+index3*3+2, nz, memory_order_relaxed);
}

inline void add_tetrahedron_2(int tetrahedron[4][3],
                              device int *level_set,
                              constant int *tetrahedron_to_triangles,
                              device int *vertex_to_index,
                              device int *vertices,
                              device atomic_int *normals,
                              device int *indices,
                              thread int *index_idx, int c[3], device int *debug_arr){
  int a1 = cube_idx_2(tetrahedron[0]);
  int a2 = cube_idx_2(tetrahedron[1]);
  int a3 = cube_idx_2(tetrahedron[2]);
  int a4 = cube_idx_2(tetrahedron[3]);
  bool b1 = level_set[a1]>80000;
  bool b2 = level_set[a2]>80000;
  bool b3 = level_set[a3]>80000;
  bool b4 = level_set[a4]>80000;
  int tetrahedron_idx = ((int)b1) + ((int)b2<<1) + ((int)b3<<2) + ((int)b4<<3);
  for (int triangle_i=0; triangle_i<12 && tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i]!=-1; triangle_i+=6){
    int triangle[3][3];
    for (int t = 0; t < 3; t++){
      char p0 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2];
      char p1 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2+1];
      int x = 0; int y = 0; int z = 0;
      if (p0==0 || p1==0){ x+=tetrahedron[0][0]; y+=tetrahedron[0][1]; z+=tetrahedron[0][2]; }
      if (p0==1 || p1==1){ x+=tetrahedron[1][0]; y+=tetrahedron[1][1]; z+=tetrahedron[1][2]; }
      if (p0==2 || p1==2){ x+=tetrahedron[2][0]; y+=tetrahedron[2][1]; z+=tetrahedron[2][2]; }
      if (p0==3 || p1==3){ x+=tetrahedron[3][0]; y+=tetrahedron[3][1]; z+=tetrahedron[3][2]; }
      triangle[t][0]=x; triangle[t][1]=y; triangle[t][2]=z;
    }
    add_triangle_2(triangle, vertex_to_index, vertices, normals, indices, index_idx,c,debug_arr);
  }
}

kernel void add_cube_2(uint3 thread_position_in_grid [[thread_position_in_grid]],
                       device int* level_set[[buffer(0)]],
                       device int* vertex_to_index[[buffer(1)]],
                       device int* indices_scan [[buffer(2)]],
                       constant int* tetrahedron_to_triangles[[buffer(3)]],
                       constant int* cube_to_tetrahedron[[buffer(4)]],
                       device int* vertices [[buffer(5)]],
                       device atomic_int* normals [[buffer(6)]],
                       device int* indices [[buffer(7)]],
                       device int* debug_arr [[buffer(8)]]){
  int c[3];
  c[0] = thread_position_in_grid.x;
  c[1] = thread_position_in_grid.y;
  c[2] = thread_position_in_grid.z;
  
  int index_idx = get_indices_idx(c,indices_scan);
  for (int tetrahedron_i = 0; tetrahedron_i<5; tetrahedron_i++){
    int tetrahedron[4][3];
    for (int i=0; i<4; i++){
      for (int j=0; j<3; j++){
        tetrahedron[i][j]=cube_to_tetrahedron[cube_to_tetrahedron_idx((c[0]+c[1]+c[2])%2,tetrahedron_i,i,j)];
      }
    }
    tetrahedron[0][0]+=c[0]; tetrahedron[0][1]+=c[1]; tetrahedron[0][2]+=c[2];
    tetrahedron[1][0]+=c[0]; tetrahedron[1][1]+=c[1]; tetrahedron[1][2]+=c[2];
    tetrahedron[2][0]+=c[0]; tetrahedron[2][1]+=c[1]; tetrahedron[2][2]+=c[2];
    tetrahedron[3][0]+=c[0]; tetrahedron[3][1]+=c[1]; tetrahedron[3][2]+=c[2];
    add_tetrahedron_2(tetrahedron, level_set, tetrahedron_to_triangles, vertex_to_index,
                      vertices, normals, indices, &index_idx, c, debug_arr);
  }
}


//Lighting:
struct Vertex {
  float4 position [[position]];
  float3 worldNormal;
  float3 worldPosition;
};

struct Uniforms {
    float4x4 modelMatrix;
    float4x4 viewProjectionMatrix;
    float3x3 normalMatrix;
};

vertex Vertex vertexShader(const device int *vertexArray [[buffer(0)]],
                           const device int *normalsArray [[buffer(1)]],
                           constant Uniforms &uniforms [[buffer(2)]],
                           unsigned int vid [[vertex_id]]){
  Vertex vertexIn;
  vertexIn.position = float4((float)vertexArray[vid*3]/N/split/2-0.5,
                             (float)vertexArray[vid*3+1]/N/split/2-0.5,
                             (float)vertexArray[vid*3+2]/N/split/2-0.5, 1);
  int a = normalsArray[vid*3];
  int b = normalsArray[vid*3+1];
  int c = normalsArray[vid*3+2];
  vertexIn.worldNormal = normalize(float3((float)normalsArray[vid*3],
                                          (float)normalsArray[vid*3+1],
                                          (float)normalsArray[vid*3+2]));
  
  float4 worldPosition = uniforms.modelMatrix * vertexIn.position;
  Vertex vertexOut;
  vertexOut.position = uniforms.viewProjectionMatrix * worldPosition;
  vertexOut.worldPosition = worldPosition.xyz;
  vertexOut.worldNormal = uniforms.normalMatrix * vertexIn.worldNormal;
  return vertexOut;
}

constant float3 ambientIntensity = 0.3;
constant float3 lightPosition(2, 2, 2); // Light position in world space
constant float3 lightColor(1, 1, 1);
constant float3 baseColor(0, 0, 1.0);

fragment float4 fragmentShader(Vertex fragmentIn [[stage_in]]){
  float3 N2 = normalize(fragmentIn.worldNormal.xyz);
  float3 L = normalize(lightPosition - fragmentIn.worldPosition.xyz);
  float3 diffuseIntensity = saturate(dot(N2, L));
  float3 finalColor = saturate(ambientIntensity + diffuseIntensity) * lightColor * baseColor;
  return float4(finalColor, 1);
}
