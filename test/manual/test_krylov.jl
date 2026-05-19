using Krylov
using Metal

A = rand(5,5)
b = rand(5)

workspace = krylov_workspace(Val(:gmres), A, b)

krylov_solve!(workspace, A, b)

# A_gpu = Metal.MtlArray(Float32.(A))
# b_gpu = Metal.MtlArray(Float32.(b))

# ws_gpu = krylov_workspace(Val(:gmres), A_gpu, b_gpu)
# krylov_solve!(ws_gpu, A_gpu, b_gpu)
# x, stats = gmres(A_gpu, b_gpu)

# appears to perform the solve on the GPU
