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
  var stopFrame: Int = 30*15
  var frameI: Int = 0
  var time: Float = 0
  var commandQueue: MTLCommandQueue!
  var pipelineState: MTLRenderPipelineState!
  var levelSet: MTLFunction!
  var levelSetPipelineState: MTLComputePipelineState!
  var marchingTetrahedra1: MTLFunction!
  var marchingTetrahedra2: MTLFunction!
  var marchingTetrahedra1PipelineState: MTLComputePipelineState!
  var marchingTetrahedra2PipelineState: MTLComputePipelineState!
  var globalScan: MTLFunction!
  var globalScanPipelineState: MTLComputePipelineState!
  var vertexBuffer: MTLBuffer!
  var normalsBuffer: MTLBuffer!
  var indicesBuffer: MTLBuffer!
  var level_set: MTLBuffer!
  var vertex_to_index: MTLBuffer!
  var pow_lookup: MTLBuffer!
  var particles_b: MTLBuffer!
  var sphere_b: MTLBuffer!
  var tetrahedron_to_triangles: MTLBuffer!
  var cube_to_tetrahedron: MTLBuffer!
  var debug_arr: MTLBuffer!
  var indices_scan: MTLBuffer!
  var device: MTLDevice!
  var depthStencilState: MTLDepthStencilState!
  var metalVideoRecorder: MetalVideoRecorder!
    
  override func viewDidLoad() {
    super.viewDidLoad()

    let mtkView = self.view as! MTKView
    device = MTLCreateSystemDefaultDevice()!
    mtkView.device = device
    mtkView.delegate = self
    commandQueue = device.makeCommandQueue()!
    let library = device.makeDefaultLibrary()!
    
    levelSet = library.makeFunction(name: "calc_level_set")
    marchingTetrahedra1 = library.makeFunction(name: "add_cube")
    marchingTetrahedra2 = library.makeFunction(name: "add_cube_2")
    globalScan = library.makeFunction(name: "globalScan")
    do {
      marchingTetrahedra1PipelineState = try device.makeComputePipelineState(function: marchingTetrahedra1);
      marchingTetrahedra2PipelineState = try device.makeComputePipelineState(function: marchingTetrahedra2);
      levelSetPipelineState = try device.makeComputePipelineState(function: levelSet);
      globalScanPipelineState = try device.makeComputePipelineState(function: globalScan);
    }catch{
      print("Error: \(error)")
    }
    
    let pipelineDescriptor = MTLRenderPipelineDescriptor()
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
    mtkView.framebufferOnly = false;
    metalVideoRecorder = MetalVideoRecorder(outputURL: "/Users/arifiorino/Movies/out.mov",
                                            size: CGSize(width: 1600, height: 1600))
    metalVideoRecorder.startRecording()
    init_animation()
    init_render()
    
    var sphere: [Int32]! = []
    for di in -split...split{
      for dj in -split...split{
        for dk in -split...split{
          if (di*di+dj*dj+dk*dk < split*split){
            sphere.append(Int32(di))
            sphere.append(Int32(dj))
            sphere.append(Int32(dk))
          }
        }
      }
    }
    debug_arr = device.makeBuffer(length: 4*100, options: MTLResourceOptions.storageModeShared)!
    indices_scan = device.makeBuffer(length: 4*(N*split)*(N*split)*(N*split), options: MTLResourceOptions.storageModeShared)!
    
    particles_b = device.makeBuffer(length: Int(n_particles) * 3 * 4, options: MTLResourceOptions.storageModeShared)!
    sphere_b = device.makeBuffer(length: Int(sphere.count) * 3 * 4, options: MTLResourceOptions.storageModeShared)!
    sphere_b.contents().copyMemory(from: sphere, byteCount: Int(sphere.count) * 3 * 4)
    
    level_set = device.makeBuffer(length: 4*(N*split+1)*(N*split+1)*(N*split+1), options: MTLResourceOptions.storageModeShared)!
    vertex_to_index = device.makeBuffer(length: 4*(N*split*2)*(N*split*2)*(N*split*2), options: MTLResourceOptions.storageModeShared)!
    let pow_lookup_a:[Int32] = [1,2,4,8,16,32,64,128,256,512,1024,2048];
    pow_lookup = device.makeBuffer(bytes: pow_lookup_a, length: 4*pow_lookup_a.count, options: MTLResourceOptions.storageModeShared)!
    
    let tetrahedron_to_triangles_a: [Int32] =
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
    let cube_to_tetrahedron_a: [Int32] =
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
    tetrahedron_to_triangles = device.makeBuffer(bytes: tetrahedron_to_triangles_a, length: 4*192, options: [])!
    cube_to_tetrahedron = device.makeBuffer(bytes: cube_to_tetrahedron_a, length: 4*120, options: [])!
    //let indices_scan = device.makeBuffer(length: 4*(N*split)*(N*split)*(N*split), options: [])!
  }
  
  func draw(in view: MTKView) {
    animate()
    //render()
    //print("actual vertices", n_vertices)
    //print("actual indices", n_indices)
    
    particles_b.contents().copyMemory(from: particles_arr, byteCount: Int(n_particles) * 3 * 4)
    
    let start = NSDate().timeIntervalSince1970
    
    let desc = MTLCommandBufferDescriptor()
    desc.errorOptions = .encoderExecutionStatus
    var commandBuffer = commandQueue.makeCommandBuffer(descriptor: desc)!
    
    var computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(levelSetPipelineState)
    computeEncoder.setBuffer(particles_b, offset: 0, index: 0)
    computeEncoder.setBuffer(level_set, offset: 0, index: 1)
    computeEncoder.setBuffer(sphere_b, offset: 0, index: 2)
    computeEncoder.dispatchThreads(MTLSizeMake(Int(n_particles),1,1), threadsPerThreadgroup: MTLSizeMake(8,1,1))
    computeEncoder.endEncoding()
    
    computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(marchingTetrahedra1PipelineState)
    computeEncoder.setBuffer(level_set, offset: 0, index: 0)
    computeEncoder.setBuffer(vertex_to_index, offset: 0, index: 1)
    computeEncoder.setBuffer(tetrahedron_to_triangles, offset: 0, index: 2)
    computeEncoder.setBuffer(cube_to_tetrahedron, offset: 0, index: 3)
    computeEncoder.setBuffer(debug_arr, offset: 0, index: 4)
    computeEncoder.setBuffer(pow_lookup, offset: 0, index: 5)
    computeEncoder.setBuffer(indices_scan, offset: 0, index: 6)
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
    
    let end = NSDate().timeIntervalSince1970
    print(end-start)
    
    var y = unsafeBitCast(vertex_to_index.contents(), to: UnsafeMutablePointer<Int32>.self)
    var y1 = Array(UnsafeBufferPointer(start: y, count: (N*split*2)*(N*split*2)*(N*split*2))) //*2
    var z = unsafeBitCast(level_set.contents(), to: UnsafeMutablePointer<Int32>.self)
    var z1 = Array(UnsafeBufferPointer(start: z, count: (N*split+1)*(N*split+1)*(N*split+1)))
    var w = unsafeBitCast(indices_scan.contents(), to: UnsafeMutablePointer<Int32>.self)
    var w1 = Array(UnsafeBufferPointer(start: w, count: (N*split)*(N*split)*(N*split)))
    
    var n_verts = Int(y1[133955070]); //HARDCODED!!!!
    print("n_vertices",n_verts)
    var n_idx = Int(w1[16777215]);
    print("n_indices",n_idx)
    
    vertexBuffer = device.makeBuffer(length: Int(n_verts) * 3 * 4, options: MTLResourceOptions.storageModeShared)!
    normalsBuffer = device.makeBuffer(length: Int(n_verts) * 3 * 4, options: MTLResourceOptions.storageModeShared)!
    indicesBuffer = device.makeBuffer(length: Int(n_idx) * 4, options: MTLResourceOptions.storageModeShared)!
    
    commandBuffer = commandQueue.makeCommandBuffer(descriptor: desc)!
    
    computeEncoder = commandBuffer.makeComputeCommandEncoder()!
    computeEncoder.setComputePipelineState(marchingTetrahedra2PipelineState)
    computeEncoder.setBuffer(level_set, offset: 0, index: 0)
    computeEncoder.setBuffer(vertex_to_index, offset: 0, index: 1)
    computeEncoder.setBuffer(indices_scan, offset: 0, index: 2)
    computeEncoder.setBuffer(tetrahedron_to_triangles, offset: 0, index: 3)
    computeEncoder.setBuffer(cube_to_tetrahedron, offset: 0, index: 4)
    computeEncoder.setBuffer(vertexBuffer, offset: 0, index: 5)
    computeEncoder.setBuffer(normalsBuffer, offset: 0, index: 6)
    computeEncoder.setBuffer(indicesBuffer, offset: 0, index: 7)
    computeEncoder.setBuffer(debug_arr, offset: 0, index: 8)
    
    computeEncoder.dispatchThreads(MTLSizeMake(N*split, N*split, N*split), threadsPerThreadgroup: MTLSizeMake(1,1,1))
    computeEncoder.endEncoding()
    
    let renderPassDescriptor = view.currentRenderPassDescriptor!
    renderPassDescriptor.colorAttachments[0].clearColor = MTLClearColorMake(1, 1, 1, 1)
    let renderEncoder = commandBuffer.makeRenderCommandEncoder(descriptor: renderPassDescriptor)!
    renderEncoder.setDepthStencilState(depthStencilState)
    renderEncoder.setRenderPipelineState(pipelineState)
    
    time += 1 / 60.0 / 2.0
    let angle = -time
    let modelMatrix = float4x4(rotationAbout: SIMD3<Float>(0, 1, 0), by: angle) *  float4x4(scaleBy: 2)

    let viewMatrix = float4x4(translationBy: SIMD3<Float>(0, 0, -2))
    let aspectRatio = Float(view.drawableSize.width / view.drawableSize.height)
    let projectionMatrix = float4x4(perspectiveProjectionFov: Float.pi / 3, aspectRatio: aspectRatio, nearZ: 0.1, farZ: 100)
    let viewProjectionMatrix = projectionMatrix * viewMatrix
    
    var uniforms = Uniforms(modelMatrix: modelMatrix, viewProjectionMatrix: viewProjectionMatrix, normalMatrix: modelMatrix.normalMatrix)
    
    renderEncoder.setVertexBytes(&uniforms, length: MemoryLayout<Uniforms>.size, index: 2)

    renderEncoder.setVertexBuffer(vertexBuffer, offset: 0, index: 0)
    renderEncoder.setVertexBuffer(normalsBuffer, offset: 0, index: 1)//Int(n_idx)
    renderEncoder.drawIndexedPrimitives(type: .triangle, indexCount: n_idx, indexType: MTLIndexType.uint32,
                                        indexBuffer: indicesBuffer, indexBufferOffset: 0)
    renderEncoder.endEncoding()
    
    if (frameI<stopFrame){
      let texture = (self.view as! MTKView).currentDrawable!.texture
      commandBuffer.addCompletedHandler { commandBuffer in
        self.metalVideoRecorder.writeFrame(forTexture: texture, frameI: self.frameI)
      }
      print("frame",frameI,"/",stopFrame)
    }else if (frameI==stopFrame){
      self.metalVideoRecorder.endRecording({});
    }
    frameI+=1;
    commandBuffer.present(view.currentDrawable!)
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    
    let v = unsafeBitCast(vertexBuffer.contents(), to: UnsafeMutablePointer<Int32>.self)
    let v1 = Array(UnsafeBufferPointer(start: v, count: n_verts*3))
    let n = unsafeBitCast(normalsBuffer.contents(), to: UnsafeMutablePointer<Int32>.self)
    let n1 = Array(UnsafeBufferPointer(start: n, count: n_verts*3))
    let i = unsafeBitCast(indicesBuffer.contents(), to: UnsafeMutablePointer<Int32>.self)
    let i1 = Array(UnsafeBufferPointer(start: i, count: n_idx))
    let x = unsafeBitCast(debug_arr.contents(), to: UnsafeMutablePointer<Int32>.self)
    let x1 = Array(UnsafeBufferPointer(start: x, count: 100))

    print("start")
    //(x1 as NSArray).write(to: URL(fileURLWithPath: "/Users/arifiorino/Downloads/debug_arr.txt"), atomically: false)
    print("end")

  }
  
  func mtkView(_ view: MTKView, drawableSizeWillChange size: CGSize) {
    
  }

}

