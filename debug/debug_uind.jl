out = zeros(3, body.ncells)
normals = pnl._calc_normals(body)
ctrl_pts = pnl._calc_controlpoints(body, normals)

pnl._rigid_wake_U!(body, ctrl_pts, out)

out2 = zeros(3, body.ncells)
pnl._rigid_wake_U2!(body, ctrl_pts, out2)

out_back = zeros(3,1)
target_back = zeros(3,1)
i = 2
pt_idx = body.cells[:,i]
target_back .= sum(body.grid._nodes[:,j] for j in pt_idx, dims=2)
target_back .+= body.Das[:,i] .* 3

pnl._rigid_wake_U2!(body, target_back, out_back)
