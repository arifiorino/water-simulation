//
//  Render.metal
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include <metal_stdlib>
using namespace metal;

struct VertexOut {
    float4 color;
    float4 pos [[position]];
    float pointsize[[point_size]];
};

vertex VertexOut vertexShader(const device float *vertexArray [[buffer(0)]], unsigned int vid [[vertex_id]]){
  VertexOut out;
  out.color = float4(0,0,1,1);
  out.pos = float4(vertexArray[vid*2]/4.0f-1.0f, vertexArray[vid*2+1]/4.0f-1.0f, 0, 1);
  out.pointsize = 5.0;
  return out;
}

fragment float4 fragmentShader(VertexOut interpolated [[stage_in]]){
  return interpolated.color;
}
