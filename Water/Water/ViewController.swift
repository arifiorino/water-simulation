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
    initAnimation()
    
  }
  
  func draw(in view: MTKView) {
    animate()
    marching_squares()
    vertexBuffer = device.makeBuffer(length: Int(n_triangles) * 3 * 2 * 4, options: [])!
    vertexBuffer.contents().copyMemory(from: triangles, byteCount: Int(n_triangles) * 3 * 2 * 4)
    let commandBuffer = commandQueue.makeCommandBuffer()!
    let renderPassDescriptor = view.currentRenderPassDescriptor!
    renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(1, 1, 1, 1)
    let renderEncoder = commandBuffer.makeRenderCommandEncoder(descriptor: renderPassDescriptor)!
    renderEncoder.setRenderPipelineState(pipelineState)
    renderEncoder.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
    renderEncoder.drawPrimitives(type: .triangle, vertexStart: 0, vertexCount: Int(n_triangles)*3)
    renderEncoder.endEncoding()
    commandBuffer.present(view.currentDrawable!)
    commandBuffer.commit()
  }
  
  func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    
  }

}

