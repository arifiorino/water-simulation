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
#define metaballsCutoff 80000

void double_scan(threadgroup int* A, threadgroup int* B, int j, constant int *powLookup){
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

int vertex_idx(int3 v){
  return  v.x + v.y*(N*split*2) + v.z*(N*split*2)*(N*split*2);
}
int level_set_idx(int3 c){
  return  c.x + c.y*(N*split+1) + c.z*(N*split+1)*(N*split+1);
}
int cube_idx(int3 c){
  return c.x + c.y*(N*split) + c.z*(N*split)*(N*split);
}
int cube_to_tetrahedron_idx(int a, int b, int c, int d){
  return d + c*3 + b*3*4 + a*3*4*5;
}

int metaballs(float3 a, float3 b){
  float s2 = distance_squared(a, b);
  if (s2 < 1)
    return (int)((1-s2)*(1-s2)*(1-s2)*100000);
  return 0;
}

kernel void calc_level_set(uint3 pos [[thread_position_in_grid]],
                           device atomic_int* level_set [[buffer(0)]],
                           device float* particles [[buffer(1)]],
                           constant int* sphere [[buffer(2)]]){
  int idx = pos.x;
  float3 particle = float3(particles[idx*3],particles[idx*3+1],particles[idx*3+2]);
  int3 center = int3(particle*split);
  for (int sphereI=0; sphereI<sphereN*3; sphereI+=3){
    int3 diff = int3(sphere[sphereI], sphere[sphereI+1], sphere[sphereI+2]);
    int3 box = center+diff;
    if (all(box>=0) && all(box < N*split+1)){
      int x = level_set_idx(box);
      int z =atomic_load_explicit(level_set+x, memory_order_relaxed);
      if (z>80000)
        continue;
      int m = metaballs(particle, float3(box)/split);
      atomic_fetch_add_explicit(level_set+x, m, memory_order_relaxed);
    }
  }
}

void add_vertex_1(int3 v, int3 c, thread int *n_vertices, device int *vertex_to_index){
  if (all(c == v / 2)){
    int vi = vertex_idx(v);
    if (vertex_to_index[vi]==-1){
        vertex_to_index[vi]=*n_vertices;
        *n_vertices+=1;
    }
  }
}

void add_triangle_1(int3 triangle[3], int3 c,
                    thread int *n_vertices, thread int *n_indices,device int *vertex_to_index){
  add_vertex_1(triangle[0], c, n_vertices, vertex_to_index);
  add_vertex_1(triangle[1], c, n_vertices, vertex_to_index);
  add_vertex_1(triangle[2], c, n_vertices, vertex_to_index);
  *n_indices+=3;
}

void add_tetrahedron_1(int3 tetrahedron[4], int3 c,
                       thread int *n_vertices, thread int *n_indices,
                       device int *level_set, device int *vertex_to_index,
                       constant int *tetrahedron_to_triangles){
  int4 a = int4(level_set[level_set_idx(tetrahedron[0])], level_set[level_set_idx(tetrahedron[1])],
                level_set[level_set_idx(tetrahedron[2])], level_set[level_set_idx(tetrahedron[3])]);
  bool4 b = a > metaballsCutoff;
  int tetrahedron_idx = ((int)b.x) + ((int)b.y<<1) + ((int)b.z<<2) + ((int)b.w<<3);
  for (int triangle_i=0; triangle_i<12 && tetrahedron_to_triangles[tetrahedron_idx*12 + triangle_i]!=-1; triangle_i+=6){
    int3 triangle[3];
    for (int t = 0; t < 3; t++){
      char p0 = tetrahedron_to_triangles[tetrahedron_idx*12 + triangle_i + t*2];
      char p1 = tetrahedron_to_triangles[tetrahedron_idx*12 + triangle_i + t*2+1];
      triangle[t] = tetrahedron[p0] + tetrahedron[p1];
    }
    add_triangle_1(triangle, c, n_vertices, n_indices, vertex_to_index);
  }
}

//Threadgroup: 32kb
kernel void add_cube_1(uint3 pos [[thread_position_in_grid]],
                       threadgroup int* tg_A [[threadgroup(0)]],
                       threadgroup int* tg_B [[threadgroup(1)]],
                       device int* level_set [[buffer(0)]],
                       device int* vertex_to_index[[buffer(1)]],
                       device int* indices_scan[[buffer(2)]],
                       constant int* cube_to_tetrahedron[[buffer(3)]],
                       constant int* tetrahedron_to_triangles[[buffer(4)]],
                       constant int* pow_lookup[[buffer(5)]]){
  int3 c = (int3)pos;
  int n_vertices=0;
  int n_indices=0;
  for (int di=0; di<2; di++)
    for (int dj=0; dj<2; dj++)
      for (int dk=0; dk<2; dk++)
        vertex_to_index[vertex_idx(c * 2 + int3(di, dj, dk))] = -1;
  
  for (int tetrahedron_i = 0; tetrahedron_i<5; tetrahedron_i++){
    int3 tetrahedron[4];
    for (int i=0; i<4; i++){
      int idx = cube_to_tetrahedron_idx((c[0]+c[1]+c[2])%2,tetrahedron_i,i,0);
      tetrahedron[i] = c + int3(cube_to_tetrahedron[idx], cube_to_tetrahedron[idx+1], cube_to_tetrahedron[idx+2]);
    }
    add_tetrahedron_1(tetrahedron, c, &n_vertices, &n_indices,
                      level_set, vertex_to_index, tetrahedron_to_triangles);
  }
  int3 g = c % 8;
  int idx = g.x + g.y*8 + g.z*64;
  tg_A[idx]=n_vertices;
  tg_B[idx]=n_indices;
  double_scan(tg_A, tg_B, idx, pow_lookup);
  vertex_to_index[vertex_idx(c * 2)]=tg_A[idx];
  indices_scan[cube_idx(c)] = tg_B[idx];
}


kernel void globalScan(device int* vertex_to_index[[buffer(0)]],
                       device int* indices_scan[[buffer(1)]]){
  int n_verts = 0;
  for (int k=14; k<N*split*2; k+=16){
    for (int j=14; j<N*split*2; j+=16){
      for (int i=14; i<N*split*2; i+=16){
        int3 v = int3(i, j, k);
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
        int3 c = int3(i, j, k);
        int idx = cube_idx(c);
        n_idx += indices_scan[idx];
        indices_scan[idx] = n_idx;
      }
    }
  }
}

//PASS 2!!!!!!!!!

int get_vertex_idx(int3 v, device int* vertex_to_index){
  int3 c = v / 2;
  int3 box = c / 8;
  int prevBox = 0;
  if (any(box != 0)){
    int idx =  box.x + box.y*32 + box.z*32*32 -  1;
    int3 pos = int3(idx%32, (idx/32)%32, idx/32/32);
    pos = 2 * (pos * 8 + 7);
    prevBox = vertex_to_index[vertex_idx(pos)];
  }
  int3 grid = c % 8;
  int prevGrid = 0;
  if (any(grid != 0)){
    int idx = grid.x + grid.y*8 + grid.z*64 - 1;
    int3 pos = int3(idx % 8, (idx/8)%8, idx/64);
    pos = 2 * (box*8 + pos);
    prevGrid = vertex_to_index[vertex_idx(pos)];
  }
  int vi = vertex_to_index[vertex_idx(v)];
  return prevBox + prevGrid + vi;
}

int get_indices_idx(int3 c, device int* indices_scan){
  int3 box = c / 8;
  int prevBox = 0;
  if (any(box != 0)){
    int idx =  box.x + box.y*32 + box.z*32*32 -  1;
    int3 pos = int3(idx%32, (idx/32)%32, idx/32/32);
    pos = pos * 8 + 7;
    prevBox = indices_scan[cube_idx(pos)];
  }
  int3 grid = c % 8;
  int prevGrid = 0;
  if (any(grid != 0)){
    int idx = grid.x + grid.y*8 + grid.z*64 - 1;
    int3 pos = int3(idx % 8, (idx/8)%8, idx/64);
    pos = box*8 + pos;
    prevGrid = indices_scan[cube_idx(pos)];
  }
  return prevBox + prevGrid;
}

int3 cross_prod(int3 a, int3 b){
  return int3(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

int add_vertex_2(int3 v, int3 c,
                 device int* vertices, device int* vertex_to_index){
  int idx = get_vertex_idx(v, vertex_to_index);
  if (all(c == v / 2)){ //It is your vertex
    vertices[idx*3  ]=v.x;
    vertices[idx*3+1]=v.y;
    vertices[idx*3+2]=v.z;
  }
  return idx;
}

void add_triangle_2(int3 triangle[3], int3 c,
                    thread int *index_idx,
                    device int *vertex_to_index,
                    device int *vertices,
                    device atomic_int *normals,
                    device int *indices){
  int3 normal = cross_prod(triangle[1]-triangle[0], triangle[2]-triangle[1]);
  for (int i=0; i<3; i++){
    int index = add_vertex_2(triangle[i], c, vertices, vertex_to_index);
    indices[*index_idx + i] = index;
    atomic_fetch_add_explicit(normals+index*3  , normal.x, memory_order_relaxed);
    atomic_fetch_add_explicit(normals+index*3+1, normal.y, memory_order_relaxed);
    atomic_fetch_add_explicit(normals+index*3+2, normal.z, memory_order_relaxed);
  }
  *index_idx+=3;
}

void add_tetrahedron_2(int3 tetrahedron[4], int3 c,
                       thread int *index_idx,
                       device int *level_set,
                       device int *vertex_to_index,
                       device int *vertices,
                       device atomic_int *normals,
                       device int *indices,
                       constant int *tetrahedron_to_triangles){
  int4 a = int4(level_set[level_set_idx(tetrahedron[0])], level_set[level_set_idx(tetrahedron[1])],
                level_set[level_set_idx(tetrahedron[2])], level_set[level_set_idx(tetrahedron[3])]);
  bool4 b = a > metaballsCutoff;
  int tetrahedron_idx = ((int)b.x) + ((int)b.y<<1) + ((int)b.z<<2) + ((int)b.w<<3);
  for (int triangle_i=0; triangle_i<12 && tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i]!=-1; triangle_i+=6){
    int3 triangle[3];
    for (int t = 0; t < 3; t++){
      char p0 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2];
      char p1 = tetrahedron_to_triangles[tetrahedron_idx*12+triangle_i + t*2+1];
      triangle[t] = tetrahedron[p0] + tetrahedron[p1];
    }
    add_triangle_2(triangle, c, index_idx, vertex_to_index, vertices, normals, indices);
  }
}

kernel void add_cube_2(uint3 thread_position_in_grid [[thread_position_in_grid]],
                       device int* vertices [[buffer(0)]],
                       device atomic_int* normals [[buffer(1)]],
                       device int* indices [[buffer(2)]],
                       device int* level_set[[buffer(3)]],
                       device int* vertex_to_index[[buffer(4)]],
                       device int* indices_scan [[buffer(5)]],
                       constant int* cube_to_tetrahedron[[buffer(6)]],
                       constant int* tetrahedron_to_triangles[[buffer(7)]]){
  int3 c = (int3)thread_position_in_grid;
  int index_idx = get_indices_idx(c,indices_scan);
  for (int tetrahedron_i = 0; tetrahedron_i<5; tetrahedron_i++){
    int3 tetrahedron[4];
    for (int i=0; i<4; i++){
      int idx = cube_to_tetrahedron_idx((c[0]+c[1]+c[2])%2,tetrahedron_i,i,0);
      tetrahedron[i].x=cube_to_tetrahedron[idx];
      tetrahedron[i].y=cube_to_tetrahedron[idx+1];
      tetrahedron[i].z=cube_to_tetrahedron[idx+2];
      tetrahedron[i]+=c;
    }
    add_tetrahedron_2(tetrahedron, c, &index_idx, level_set, vertex_to_index,
                      vertices, normals, indices, tetrahedron_to_triangles);
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

vertex Vertex vertexShader(unsigned int vid [[vertex_id]],
                           const device int *vertexArray [[buffer(0)]],
                           const device int *normalsArray [[buffer(1)]],
                           constant Uniforms &uniforms [[buffer(2)]]){
  Vertex vertexIn;
  float3 pos = float3(vertexArray[vid*3], vertexArray[vid*3+1], vertexArray[vid*3+2]);
  pos = pos/N/split/2 - 0.5;
  vertexIn.position = float4(pos, 1);
  vertexIn.worldNormal = normalize(float3(normalsArray[vid*3], normalsArray[vid*3+1], normalsArray[vid*3+2]));
  float4 worldPosition = uniforms.modelMatrix * vertexIn.position;
  Vertex vertexOut;
  vertexOut.position = uniforms.viewProjectionMatrix * worldPosition;
  vertexOut.worldPosition = worldPosition.xyz;
  vertexOut.worldNormal = uniforms.normalMatrix * vertexIn.worldNormal;
  return vertexOut;
}

constant float3 ambientIntensity = 0.3;
constant float3 lightPosition(2, 2, 2);
constant float3 lightColor(1, 1, 1);
constant float3 baseColor(0, 0, 1.0);

fragment float4 fragmentShader(Vertex fragmentIn [[stage_in]]){
  float3 N2 = normalize(fragmentIn.worldNormal.xyz);
  float3 L = normalize(lightPosition - fragmentIn.worldPosition.xyz);
  float3 diffuseIntensity = saturate(dot(N2, L));
  float3 finalColor = saturate(ambientIntensity + diffuseIntensity) * lightColor * baseColor;
  return float4(finalColor, 1);
}
