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
  float3 normal;
};

vertex Vertex vertexShader(const device float *vertexArray [[buffer(0)]],
                           const device float *normalsArray [[buffer(1)]],
                           const device float4x4 *projMatrix [[buffer(2)]],
                           unsigned int vid [[vertex_id]]){
  Vertex out;
  out.pos = float4(vertexArray[vid*3], vertexArray[vid*3+1], vertexArray[vid*3+2], 1);
  out.pos = (*projMatrix) * out.pos;
  out.normal = float3(normalsArray[vid*3], normalsArray[vid*3+1], normalsArray[vid*3+2]);
  return out;
}

constant float3 ambientIntensity = 0.4;
constant float3 lightPosition(0, -2, 4); // Light position in world space
constant float3 lightColor(1, 1, 1);
constant float3 worldCameraPosition(0, 0, 0);
constant float3 baseColor(0, 0, 1);
constant float specularPower = 200;

fragment float4 fragmentShader(Vertex fragmentIn [[stage_in]]){
  float3 N = normalize(fragmentIn.normal);
  float3 pos = float3(fragmentIn.pos.x, fragmentIn.pos.y, fragmentIn.pos.z);
  float3 L = normalize(lightPosition - pos);
  float3 diffuseIntensity = saturate(dot(N, L));
  float3 V = normalize(worldCameraPosition - pos);
  float3 H = normalize(L + V);
  float specularBase = saturate(dot(N, H));
  float specularIntensity = powr(specularBase, specularPower);
  float3 finalColor = saturate(ambientIntensity + diffuseIntensity) * baseColor * lightColor +
                      specularIntensity * lightColor;
  return float4(finalColor, 1);
}
