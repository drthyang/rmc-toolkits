// WebGPU-accelerated Gaussian KDE density map.
//
// The density grid is an embarrassingly parallel map: every cell density[y][x]
// is an independent sum of an anisotropic 2D Gaussian over the sample set. This
// module runs that map as a compute shader (one invocation per cell) inside the
// KDE worker, where navigator.gpu is available on dedicated workers.
//
// It is a drop-in accelerator, not a replacement: computeDensityGpu returns the
// same grid x grid array the CPU loop produces, or null/throws when WebGPU is
// unavailable or fails. The caller (localKdeWorker.js) falls back to the JS loop
// in those cases, so devices without WebGPU behave exactly as before.

const KDE_WGSL = `
struct Params {
    inv    : vec4<f32>,   // inv00, inv01, inv11, normalizer
    origin : vec4<f32>,   // xMin, yMin, xStep, yStep
    dims   : vec4<u32>,   // grid, sampleCount, pad, pad
};

@group(0) @binding(0) var<uniform> P : Params;
@group(0) @binding(1) var<storage, read>       samples : array<vec2<f32>>;
@group(0) @binding(2) var<storage, read_write> density : array<f32>;

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) gid : vec3<u32>) {
    let grid = P.dims.x;
    if (gid.x >= grid || gid.y >= grid) {
        return;
    }
    let gx = P.origin.x + f32(gid.x) * P.origin.z;
    let gy = P.origin.y + f32(gid.y) * P.origin.w;
    let n = P.dims.y;
    var sum : f32 = 0.0;
    for (var i : u32 = 0u; i < n; i = i + 1u) {
        let s = samples[i];
        let dx = gx - s.x;
        let dy = gy - s.y;
        let e = -0.5 * (P.inv.x * dx * dx + 2.0 * P.inv.y * dx * dy + P.inv.z * dy * dy);
        if (e > -60.0) {
            sum = sum + exp(e);
        }
    }
    density[gid.y * grid + gid.x] = sum * P.inv.w;
}
`;

// GPU setup, buffer uploads, and the mapAsync readback round-trip carry a fixed
// per-call cost; below this many work units (cells x samples) the JS loop wins.
const GPU_MIN_WORK = 2_000_000;

export const shouldUseGpu = (grid, sampleCount) => grid * grid * sampleCount >= GPU_MIN_WORK;

// Acquire the adapter/device once per worker lifetime and cache the resolved
// promise. Any failure (no navigator.gpu, no adapter, requestDevice or pipeline
// rejection) collapses to null so callers fall back to the CPU loop and we never
// retry init on every message. A lost device clears the cache so a later message
// can re-initialize.
let gpuInitPromise = null;

const getGpu = () => {
    if (gpuInitPromise) return gpuInitPromise;
    gpuInitPromise = (async () => {
        if (typeof navigator === 'undefined' || !navigator.gpu) return null;
        const adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
        if (!adapter) return null;
        const device = await adapter.requestDevice();
        const module = device.createShaderModule({ code: KDE_WGSL });
        const pipeline = device.createComputePipeline({
            layout: 'auto',
            compute: { module, entryPoint: 'main' }
        });
        device.lost.then(() => {
            gpuInitPromise = null;
        });
        return { device, pipeline };
    })().catch(() => null);
    return gpuInitPromise;
};

// Compute the density grid on the GPU. Resolves to a grid x grid array of plain
// JS numbers (matching the CPU path) or null when WebGPU is unavailable; throws
// only on an unexpected GPU runtime error, which the caller treats as a fallback.
export const computeDensityGpu = async ({ samples, kernel, grid, xMin, yMin, xStep, yStep }) => {
    const gpu = await getGpu();
    if (!gpu) return null;
    const { device, pipeline } = gpu;
    const sampleCount = samples.length;

    // Samples: tightly-packed [x0, y0, x1, y1, ...] -> array<vec2<f32>> (8-byte stride).
    const sampleData = new Float32Array(sampleCount * 2);
    for (let i = 0; i < sampleCount; i += 1) {
        sampleData[2 * i] = samples[i][0];
        sampleData[2 * i + 1] = samples[i][1];
    }
    const sampleBuffer = device.createBuffer({
        size: Math.max(sampleData.byteLength, 16),
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST
    });
    device.queue.writeBuffer(sampleBuffer, 0, sampleData);

    // Params: three vec4 lanes (48 bytes). Packing into vec4s sidesteps uniform
    // alignment rules; the dims lane is read through a Uint32 view of the same buffer.
    const paramData = new ArrayBuffer(48);
    const paramFloats = new Float32Array(paramData);
    const paramUints = new Uint32Array(paramData);
    paramFloats[0] = kernel.inv00;
    paramFloats[1] = kernel.inv01;
    paramFloats[2] = kernel.inv11;
    paramFloats[3] = kernel.normalizer;
    paramFloats[4] = xMin;
    paramFloats[5] = yMin;
    paramFloats[6] = xStep;
    paramFloats[7] = yStep;
    paramUints[8] = grid;
    paramUints[9] = sampleCount;
    const paramBuffer = device.createBuffer({
        size: paramData.byteLength,
        usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
    });
    device.queue.writeBuffer(paramBuffer, 0, paramData);

    // Output storage buffer + a mappable readback buffer.
    const outputSize = grid * grid * 4;
    const outputBuffer = device.createBuffer({
        size: outputSize,
        usage: GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_SRC
    });
    const readBuffer = device.createBuffer({
        size: outputSize,
        usage: GPUBufferUsage.MAP_READ | GPUBufferUsage.COPY_DST
    });

    const bindGroup = device.createBindGroup({
        layout: pipeline.getBindGroupLayout(0),
        entries: [
            { binding: 0, resource: { buffer: paramBuffer } },
            { binding: 1, resource: { buffer: sampleBuffer } },
            { binding: 2, resource: { buffer: outputBuffer } }
        ]
    });

    const encoder = device.createCommandEncoder();
    const pass = encoder.beginComputePass();
    pass.setPipeline(pipeline);
    pass.setBindGroup(0, bindGroup);
    const workgroups = Math.ceil(grid / 8);
    pass.dispatchWorkgroups(workgroups, workgroups);
    pass.end();
    encoder.copyBufferToBuffer(outputBuffer, 0, readBuffer, 0, outputSize);
    device.queue.submit([encoder.finish()]);

    try {
        await readBuffer.mapAsync(GPUMapMode.READ);
        const flat = new Float32Array(readBuffer.getMappedRange()).slice();
        readBuffer.unmap();

        // Reshape the flat buffer into density[y][x] of plain numbers so the rest
        // of the worker (log scaling, min/max, contours) is identical to the CPU path.
        const density = new Array(grid);
        for (let y = 0; y < grid; y += 1) {
            const row = new Array(grid);
            const base = y * grid;
            for (let x = 0; x < grid; x += 1) {
                row[x] = flat[base + x];
            }
            density[y] = row;
        }
        return density;
    } finally {
        sampleBuffer.destroy();
        paramBuffer.destroy();
        outputBuffer.destroy();
        readBuffer.destroy();
    }
};
