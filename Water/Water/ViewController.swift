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

let N = 32;
let split = 8;

class ViewController: NSViewController, MTKViewDelegate {
  var stopFrame: Int = 30*20
  var frameI: Int = 0
  var angle: Float = 0
  
  var device: MTLDevice!
  var commandQueue: MTLCommandQueue!
  var levelSetPipelineState: MTLComputePipelineState!
  var march1PipelineState: MTLComputePipelineState!
  var globalScanPipelineState: MTLComputePipelineState!
  var march2PipelineState: MTLComputePipelineState!
  var depthStencilState: MTLDepthStencilState!
  var renderPipelineState: MTLRenderPipelineState!
  
  var vertices: MTLBuffer!
  var normals: MTLBuffer!
  var indices: MTLBuffer!
  var level_set: MTLBuffer!
  var vertex_to_index: MTLBuffer!
  var pow_lookup: MTLBuffer!
  var particles: MTLBuffer!
  var sphere: MTLBuffer!
  var tetrahedron_to_triangles: MTLBuffer!
  var cube_to_tetrahedron: MTLBuffer!
  var indices_scan: MTLBuffer!

  var metalVideoRecorder: MetalVideoRecorder!
    
  override func viewDidLoad() {
    super.viewDidLoad()

    let mtkView = self.view as! MTKView
    device = MTLCreateSystemDefaultDevice()!
    mtkView.device = device
    mtkView.delegate = self
    commandQueue = device.makeCommandQueue()!
    let library = device.makeDefaultLibrary()!
    
    let levelSet = library.makeFunction(name: "calc_level_set")
    let march1 = library.makeFunction(name: "add_cube_1")
    let globalScan = library.makeFunction(name: "globalScan")
    let march2 = library.makeFunction(name: "add_cube_2")
    do {
      levelSetPipelineState = try device.makeComputePipelineState(function: levelSet!);
      march1PipelineState = try device.makeComputePipelineState(function: march1!);
      globalScanPipelineState = try device.makeComputePipelineState(function: globalScan!);
      march2PipelineState = try device.makeComputePipelineState(function: march2!);
    }catch{
      print("Error: \(error)")
    }
    let render = MTLRenderPipelineDescriptor()
    render.vertexFunction = library.makeFunction(name: "vertexShader")
    render.fragmentFunction = library.makeFunction(name: "fragmentShader")
    render.colorAttachments[0].pixelFormat = mtkView.colorPixelFormat
    render.depthAttachmentPixelFormat = .depth32Float
    do {
      renderPipelineState = try device.makeRenderPipelineState(descriptor: render)
    }catch{
      print("Error: \(error)")
    }
    
    let depth = MTLDepthStencilDescriptor()
    depth.depthCompareFunction = .less
    depth.isDepthWriteEnabled = true
    depthStencilState = device.makeDepthStencilState(descriptor: depth)!
    mtkView.colorPixelFormat = .bgra8Unorm
    mtkView.depthStencilPixelFormat = .depth32Float
    mtkView.framebufferOnly = false;
    
    init_animation()
    
    particles = device.makeBuffer(length: 3 * 4 * Int(n_particles), options: [])!
    var sphere_arr: [Int32]! = []
    for dj in -split...split{
      for dk in -split...split{
        for di in -split...split{
          if (di*di+dj*dj+dk*dk < split*split){
            sphere_arr.append(Int32(di))
            sphere_arr.append(Int32(dj))
            sphere_arr.append(Int32(dk))
          }
        }
      }
    }
    sphere = device.makeBuffer(bytes: sphere_arr, length: 3 * 4 * sphere_arr.count, options: [])!
    vertex_to_index = device.makeBuffer(length: 4*(N*split*2)*(N*split*2)*(N*split*2), options: [])!
    indices_scan = device.makeBuffer(length: 4*(N*split)*(N*split)*(N*split), options: [])!
    let pow_lookup_arr:[Int32] = [1,2,4,8,16,32,64,128,256,512,1024,2048];
    pow_lookup = device.makeBuffer(bytes: pow_lookup_arr, length: 4 * pow_lookup_arr.count, options: [])!
    let tetrahedron_to_triangles_arr: [Int32] =
    [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
      0, 1, 0, 3, 0, 2,-1,-1,-1,-1,-1,-1,
      0, 1, 1, 2, 1, 3,-1,-1,-1,-1,-1,-1,
      0, 2, 1, 3, 0, 3, 0, 2, 1, 2, 1, 3,
      0, 2, 2, 3, 1, 2,-1,-1,-1,-1,-1,-1,
      0, 1, 0, 3, 2, 3, 0, 1, 2, 3, 1, 2,
      0, 1, 0, 2, 2, 3, 0, 1, 2, 3, 1, 3,
      0, 3, 2, 3, 1, 3,-1,-1,-1,-1,-1,-1,
      0, 3, 1, 3, 2, 3,-1,-1,-1,-1,-1,-1,
      0, 1, 2, 3, 0, 2, 0, 1, 1, 3, 2, 3,
      0, 1, 2, 3, 0, 3, 0, 1, 1, 2, 2, 3,
      0, 2, 1, 2, 2, 3,-1,-1,-1,-1,-1,-1,
      0, 2, 0, 3, 1, 3, 0, 2, 1, 3, 1, 2,
      0, 1, 1, 3, 1, 2,-1,-1,-1,-1,-1,-1,
      0, 1, 0, 2, 0, 3,-1,-1,-1,-1,-1,-1,
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
    tetrahedron_to_triangles = device.makeBuffer(bytes: tetrahedron_to_triangles_arr,
                                                 length: 4*tetrahedron_to_triangles_arr.count,
                                                 options: [])!
    let cube_to_tetrahedron_arr: [Int32] =
     [0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0,
      1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0,
      1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0,
      1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0,
      0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1,
      1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
      1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0,
      1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0];
    cube_to_tetrahedron = device.makeBuffer(bytes: cube_to_tetrahedron_arr,
                                            length: 4*cube_to_tetrahedron_arr.count,
                                            options: [])!
    
    metalVideoRecorder = MetalVideoRecorder(outputURL: "/Users/arifiorino/Movies/out.mov",
                                            size: CGSize(width: 1600, height: 1600))
    metalVideoRecorder.startRecording()
  }
  
  func draw(in view: MTKView) {
    animate()
    
    let start = NSDate().timeIntervalSince1970
    level_set = device.makeBuffer(length: 4*(N*split+1)*(N*split+1)*(N*split+1), options: [])!
    particles.contents().copyMemory(from: particles_arr, byteCount: 3 * 4 * Int(n_particles))
    
    var commandBuffer = commandQueue.makeCommandBuffer()!
    
    var computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(levelSetPipelineState)
    computeEncoder.setBuffer(level_set, offset: 0, index: 0)
    computeEncoder.setBuffer(particles, offset: 0, index: 1)
    computeEncoder.setBuffer(sphere, offset: 0, index: 2)
    computeEncoder.dispatchThreads(MTLSizeMake(Int(n_particles),1,1), threadsPerThreadgroup: MTLSizeMake(8,1,1))
    computeEncoder.endEncoding()
    
    computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(march1PipelineState)
    computeEncoder.setBuffer(level_set, offset: 0, index: 0)
    computeEncoder.setBuffer(vertex_to_index, offset: 0, index: 1)
    computeEncoder.setBuffer(indices_scan, offset: 0, index: 2)
    computeEncoder.setBuffer(cube_to_tetrahedron, offset: 0, index: 3)
    computeEncoder.setBuffer(tetrahedron_to_triangles, offset: 0, index: 4)
    computeEncoder.setBuffer(pow_lookup, offset: 0, index: 5)
    computeEncoder.setThreadgroupMemoryLength(4*512, index: 0)
    computeEncoder.setThreadgroupMemoryLength(4*512, index: 1)
    computeEncoder.dispatchThreads(MTLSizeMake(N*split, N*split, N*split), threadsPerThreadgroup: MTLSizeMake(8,8,8))
    computeEncoder.endEncoding()
    
    computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(globalScanPipelineState)
    computeEncoder.setBuffer(vertex_to_index, offset: 0, index: 0)
    computeEncoder.setBuffer(indices_scan, offset: 0, index: 1)
    computeEncoder.dispatchThreads(MTLSizeMake(1,1,1), threadsPerThreadgroup: MTLSizeMake(1,1,1))
    computeEncoder.endEncoding()

    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    
    
    let vertex_to_index_ptr = vertex_to_index.contents().assumingMemoryBound(to: Int32.self)
    let indices_scan_ptr = indices_scan.contents().assumingMemoryBound(to: Int32.self)
    let n_verts = Int(vertex_to_index_ptr[133955070]); //HARDCODED!!!!
    let n_idx = Int(indices_scan_ptr[16777215]);
    //print("n_vertices",n_verts,"n_indices",n_idx)
    vertices = device.makeBuffer(length:  3 * 4 * Int(n_verts), options: [])!
    normals = device.makeBuffer(length: 3 * 4 * Int(n_verts), options: [])!
    indices = device.makeBuffer(length: 4 * Int(n_idx), options: [])!
    
    commandBuffer = commandQueue.makeCommandBuffer()!
    computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(march2PipelineState)
    computeEncoder.setBuffer(vertices, offset: 0, index: 0)
    computeEncoder.setBuffer(normals, offset: 0, index: 1)
    computeEncoder.setBuffer(indices, offset: 0, index: 2)
    computeEncoder.setBuffer(level_set, offset: 0, index: 3)
    computeEncoder.setBuffer(vertex_to_index, offset: 0, index: 4)
    computeEncoder.setBuffer(indices_scan, offset: 0, index: 5)
    computeEncoder.setBuffer(cube_to_tetrahedron, offset: 0, index: 6)
    computeEncoder.setBuffer(tetrahedron_to_triangles, offset: 0, index: 7)
    computeEncoder.dispatchThreads(MTLSizeMake(N*split, N*split, N*split), threadsPerThreadgroup: MTLSizeMake(1,1,1))
    computeEncoder.endEncoding()
    
    let renderPassDescriptor = view.currentRenderPassDescriptor!
    renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(1, 1, 1, 1)
    let renderEncoder = commandBuffer.makeRenderCommandEncoder(descriptor: renderPassDescriptor)!
    renderEncoder.setDepthStencilState(depthStencilState)
    renderEncoder.setRenderPipelineState(renderPipelineState)
    angle -= 1 / 60.0 / 2
    let modelMatrix = float4x4(rotationAbout: SIMD3<Float>(0, 1, 0), by: angle) *  float4x4(scaleBy: 2)
    let viewMatrix = float4x4(translationBy: SIMD3<Float>(0, 0, -2))
    let aspectRatio = Float(view.drawableSize.width / view.drawableSize.height)
    let projectionMatrix = float4x4(perspectiveProjectionFov: Float.pi / 3, aspectRatio: aspectRatio, nearZ: 0.1, farZ: 100)
    let viewProjectionMatrix = projectionMatrix * viewMatrix
    var uniforms = Uniforms(modelMatrix: modelMatrix, viewProjectionMatrix: viewProjectionMatrix, normalMatrix: modelMatrix.normalMatrix)
    
    renderEncoder.setVertexBuffer(vertices, offset: 0, index: 0)
    renderEncoder.setVertexBuffer(normals, offset: 0, index: 1)
    renderEncoder.setVertexBytes(&uniforms, length: MemoryLayout<Uniforms>.size, index: 2)
    renderEncoder.drawIndexedPrimitives(type: .triangle, indexCount: n_idx, indexType: MTLIndexType.uint32,
                                        indexBuffer: indices, indexBufferOffset: 0)
    renderEncoder.endEncoding()
    
    /*if (frameI<stopFrame){
      let texture = (self.view as! MTKView).currentDrawable!.texture
      commandBuffer.addCompletedHandler { commandBuffer in
        self.metalVideoRecorder.writeFrame(forTexture: texture, frameI: self.frameI)
      }
      print("frame",frameI,"/",stopFrame)
    }else if (frameI==stopFrame){
      self.metalVideoRecorder.endRecording({});
    }
    frameI+=1;*/
    
    commandBuffer.present(view.currentDrawable!)
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()

    let end = NSDate().timeIntervalSince1970
    print("Render: ",end-start)

  }
  
  func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    
  }

}

