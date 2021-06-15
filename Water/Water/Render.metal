//
//  Render.metal
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include <metal_stdlib>
using namespace metal;

struct Vertex {
    float4 pos [[position]];
};

vertex Vertex vertexShader(const device float *vertexArray [[buffer(0)]], unsigned int vid [[vertex_id]]){
  int N = 64;
  Vertex out;
  out.pos = float4(vertexArray[vid*3]/(N/2.0f)-1.0f,
                   vertexArray[vid*3+1]/(N/2.0f)-1.0f,
                   vertexArray[vid*3+2]/N, 1);
  return out;
}

fragment float4 fragmentShader(Vertex interpolated [[stage_in]]){
  return float4(0,0,1,1);
}
