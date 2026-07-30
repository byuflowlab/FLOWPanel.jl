lines = []
with open("examples/duct_unsteady.jl") as f:
    lines = f.readlines()

new_lines = []
for i, line in enumerate(lines):
    new_lines.append(line)
    if "xs *= d*aspectratio" in line:
        new_lines.append("""
# Remove any duplicate points to prevent zero-area panels (NaNs in FMM matrices)
idx_to_keep = [1]
for i in 2:length(xs)
    if abs(xs[i] - xs[i-1]) > 1e-12 || abs(ys[i] - ys[i-1]) > 1e-12
        push!(idx_to_keep, i)
    end
end
xs = xs[idx_to_keep]
ys = ys[idx_to_keep]
""")

with open("examples/duct_unsteady.jl", "w") as f:
    f.writelines(new_lines)
