

using Printf
using LinearAlgebra
import FLOWPanel as pnl

println("FLOWPanel loaded!")

# generate a simple body of revolution
points = [0.0 0.0; 0.5 0.5; 1.0 0.0]
shedding = zeros(Int, 6, 0)
bodytype = pnl.RigidWakeBody{pnl.VortexRing}
body = pnl.generate_revolution_liftbody(bodytype, points, 10; loop_dim=2, axis_angle=90, closed_contour=false, overwrite_shedding=shedding)
system = body

# freestream
Umag = 10.0
alpha = 0.0
Uinf(t) = Umag * [cos(alpha), 0, sin(alpha)]

# frames
frames = pnl.ReferenceFrame(system; origin=[0.0, 0.0, 0.0], v=[1.0, 0.0, 0.0]) # moving 1 m/s in x direction

# simulation
maneuver(frames, system, wake, t) = false

n_steps = 1
dt = 0.1
t_range = 0:dt:(n_steps*dt)

println("Initial node[1] x: ", body.nodes[1, 1])

# Run simulate!
wake = pnl.simulate!(system, frames, maneuver, Uinf, t_range)

println("Final node[1] x: ", body.nodes[1, 1])
