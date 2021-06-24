//
//  ViewController.swift
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

import Cocoa
import Metal
import MetalKit

struct Uniforms {
    var modelMatrix: float4x4
    var viewProjectionMatrix: float4x4
    var normalMatrix: float3x3
}

class ViewController: NSViewController, MTKViewDelegate {
  var time: Float = 0
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
    
    time += 1 / Float(view.preferredFramesPerSecond)
    let angle = -time
    let modelMatrix = float4x4(rotationAbout: SIMD3<Float>(0, 1, 0), by: angle) *  float4x4(scaleBy: 2)

    let viewMatrix = float4x4(translationBy: SIMD3<Float>(0, 0, -2))
    let aspectRatio = Float(view.drawableSize.width / view.drawableSize.height)
    let projectionMatrix = float4x4(perspectiveProjectionFov: Float.pi / 3, aspectRatio: aspectRatio, nearZ: 0.1, farZ: 100)
    let viewProjectionMatrix = projectionMatrix * viewMatrix
    
    var uniforms = Uniforms(modelMatrix: modelMatrix, viewProjectionMatrix: viewProjectionMatrix, normalMatrix: modelMatrix.normalMatrix)
    
    renderEncoder.setVertexBytes(&uniforms, length: MemoryLayout<Uniforms>.size, index: 2)

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

