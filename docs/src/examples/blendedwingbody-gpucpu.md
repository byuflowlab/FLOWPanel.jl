# GPU and CPU Acceleration

## CPU Multi-Threading

The kernels and solvers implemented in FLOWPanel are parallelized (threaded)
in CPU by default.
However, in order to activate the CPU parallelization, the user needs to
[launch Julia with multi-threading activated](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads).
For instance, to launch Julia with 4 threads:
```bash
$ julia --threads 4
```

You can then verify that the 4 threads became available:

```julia-repl
julia> Threads.nthreads()
4
```

## Porting to GPU

The solver can be seamlessly ported to GPU by indicating the type of
array to be used internally.
[The Julia GPU interface](https://juliagpu.org/) is the same for any GPU hardware and platform
(NVIDIA CUDA, AMD ROCm, and Mac Metal), however, we have only tested NVIDIA
GPUs.

For an NVIDIA GPU, first import the CUDA package before running the code of
the previous section,
```julia-repl
julia> import CUDA
```
check that the GPU hardware is ready to be used,
```julia-repl
julia> CUDA.functional()
true
```
and instead of letting the solver default its internal arrays to CPU, change
the solver call from
```julia-repl
julia> pnl.solve(body, Uinfs, Das, Dbs)
```
to
```julia-repl
julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=CUDA.CuArray{Float32})
```

For AMD GPU:
```julia-repl
julia> import AMDGPU
julia> AMDGPU.functional()
true
julia> AMDGPU.functional(:MIOpen)
true
julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=AMDGPU.ROCArray{Float32})
```

For Metal GPU:
```julia-repl
julia> import Metal
julia> Metal.functional()
true
julia> pnl.solve(body, Uinfs, Das, Dbs; GPUArray=Metal.MtlArray{Float32})
```


!!! info "GPU"
    We have only tested NVIDIA GPUs

