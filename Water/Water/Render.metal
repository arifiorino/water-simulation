//
//  Render.metal
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#include <metal_stdlib>
using namespace metal;

struct Vertex {
  float4 position [[position]];
  float3 worldNormal;
  float3 worldPosition;
  float2 texCoords;
};

struct Uniforms {
    float4x4 modelMatrix;
    float4x4 viewProjectionMatrix;
    float3x3 normalMatrix;
};

vertex Vertex vertexShader(const device float *vertexArray [[buffer(0)]],
                           const device float *normalsArray [[buffer(1)]],
                           constant Uniforms &uniforms [[buffer(2)]],
                           unsigned int vid [[vertex_id]]){
  float N=64;
  Vertex vertexIn;
  vertexIn.position = float4(vertexArray[vid*3]/N-0.5, vertexArray[vid*3+1]/N-0.5, vertexArray[vid*3+2]/N-0.5, 1);
  vertexIn.worldNormal = normalize(float3(normalsArray[vid*3], normalsArray[vid*3+1], normalsArray[vid*3+2]));
  float4 worldPosition = uniforms.modelMatrix * vertexIn.position;
  Vertex vertexOut;
  vertexOut.position = uniforms.viewProjectionMatrix * worldPosition;
  vertexOut.worldPosition = worldPosition.xyz;
  vertexOut.worldNormal = uniforms.normalMatrix * vertexIn.worldNormal;
  vertexOut.texCoords = vertexIn.texCoords;
  return vertexOut;
}

constant float3 ambientIntensity = 0.1;
constant float3 lightPosition(2, 2, 2); // Light position in world space
constant float3 lightColor(1, 1, 1);
constant float3 baseColor(1.0, 0, 0);

fragment float4 fragmentShader(Vertex fragmentIn [[stage_in]]){
  float3 N = normalize(fragmentIn.worldNormal.xyz);
  float3 L = normalize(lightPosition - fragmentIn.worldPosition.xyz);
  float3 diffuseIntensity = saturate(dot(N, L));
  float3 finalColor = saturate(ambientIntensity + diffuseIntensity) * lightColor * baseColor;
  return float4(finalColor, 1);
}
