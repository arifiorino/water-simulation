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
    do {
      pipelineState = try device.makeRenderPipelineState(descriptor: pipelineDescriptor)
    }catch{
      print("Error: \(error)")
    }
    
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
    renderEncoder.setRenderPipelineState(pipelineState)
    renderEncoder.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
    renderEncoder.setVertexBuffer(normalsBuffer, offset: 0, index: 1)
    let viewMatrix = float4x4(translationBy: simd_float3(-0.5, -0.5, -1)) * float4x4(scaleBy: 1/14);
    let projectionMatrix = float4x4(perspectiveProjectionFov: Float.pi / 3, aspectRatio: 1, nearZ: 0.1, farZ: 1)
    var viewProjectionMatrix = projectionMatrix * viewMatrix
    renderEncoder.setVertexBytes(&viewProjectionMatrix, length: MemoryLayout<float4x4>.size, index: 2)
    renderEncoder.drawIndexedPrimitives(type: .triangle, indexCount: Int(n_indices), indexType: MTLIndexType.uint32,
                                        indexBuffer: indicesBuffer, indexBufferOffset: 0)
    renderEncoder.endEncoding()
    commandBuffer.present(view.currentDrawable!)
    commandBuffer.commit()
    
  }
  
  func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    
  }

}

