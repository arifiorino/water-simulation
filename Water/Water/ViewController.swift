//
//  ViewController.swift
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

import Cocoa
import Metal
import MetalKit

class ViewController: NSViewController, MTKViewDelegate {
    
  var commandQueue: MTLCommandQueue!
  var pipelineState: MTLRenderPipelineState!
  var vertexBuffer: MTLBuffer!
  var normalsBuffer: MTLBuffer!
  var indicesBuffer: MTLBuffer!
  var device: MTLDevice!
  var depthStencilState: MTLDepthStencilState!
    
  override func viewDidLoad() {
    super.viewDidLoad()

    let mtkView = self.view as! MTKView
    device = MTLCreateSystemDefaultDevice()!
    mtkView.device = device
    mtkView.delegate = self
    commandQueue = device.makeCommandQueue()!
    let pipelineDescriptor = MTLRenderPipelineDescriptor()
    let library = device.makeDefaultLibrary()!
    pipelineDescriptor.vertexFunction = library.makeFunction(name: "vertexShader")
    pipelineDescriptor.fragmentFunction = library.makeFunction(name: "fragmentShader")
    pipelineDescriptor.colorAttachments[0].pixelFormat = mtkView.colorPixelFormat
    pipelineDescriptor.depthAttachmentPixelFormat = .depth32Float
    do {
      pipelineState = try device.makeRenderPipelineState(descriptor: pipelineDescriptor)
    }catch{
      print("Error: \(error)")
    }
    
    let depthStencilDescriptor = MTLDepthStencilDescriptor()
    depthStencilDescriptor.depthCompareFunction = .less
    depthStencilDescriptor.isDepthWriteEnabled = true
    depthStencilState = device.makeDepthStencilState(descriptor: depthStencilDescriptor)!
    mtkView.colorPixelFormat = .bgra8Unorm
    mtkView.depthStencilPixelFormat = .depth32Float
    
    
    init_animation()
  }
  
  func draw(in view: MTKView) {
    animate()
    render()
    vertexBuffer = device.makeBuffer(length: Int(n_vertices) * 3 * 4, options: [])!
    vertexBuffer.contents().copyMemory(from: vertices, byteCount: Int(n_vertices) * 3 * 4)
    normalsBuffer = device.makeBuffer(length: Int(n_vertices) * 3 * 4, options: [])!
    normalsBuffer.contents().copyMemory(from: normals, byteCount: Int(n_vertices) * 3 * 4)
    indicesBuffer = device.makeBuffer(length: Int(n_indices) * 4, options: [])!
    indicesBuffer.contents().copyMemory(from: indices, byteCount: Int(n_indices) * 4)
    let commandBuffer = commandQueue.makeCommandBuffer()!
    let renderPassDescriptor = view.currentRenderPassDescriptor!
    renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(1, 1, 1, 1)
    let renderEncoder = commandBuffer.makeRenderCommandEncoder(descriptor: renderPassDescriptor)!
    renderEncoder.setDepthStencilState(depthStencilState)
    renderEncoder.setRenderPipelineState(pipelineState)
    
    let N: Float = 32.0;
    let scaleMatrix = float4x4(scaleByX: 1.0/N*2, scaleByY: 1.0/N*2, scaleByZ: 1.0/N*2);
    let translationMatrix = float4x4(translationBy: simd_float3(-1, -1, -5));
    let projectionMatrix = float4x4(perspectiveProjectionFov: Float.pi / 6, aspectRatio: 1, nearZ: 0.1, farZ: 10)
    var viewProjectionMatrix = projectionMatrix * translationMatrix * scaleMatrix
    renderEncoder.setVertexBytes(&viewProjectionMatrix, length: MemoryLayout<float4x4>.size, index: 2)
    
    renderEncoder.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
    renderEncoder.setVertexBuffer(normalsBuffer, offset: 0, index: 1)
    renderEncoder.drawIndexedPrimitives(type: .triangle, indexCount: Int(n_indices), indexType: MTLIndexType.uint32,
                                        indexBuffer: indicesBuffer, indexBufferOffset: 0)
    renderEncoder.endEncoding()
    commandBuffer.present(view.currentDrawable!)
    commandBuffer.commit()
    
  }
  
  func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    
  }

}

