doc=
"""
usage: 
stencil.jl <outputdir> <builddir> <refine> <scheme> <vtk>

     builddir path      Directory to write discretization.
     outputdir path     Directory to write output data to.
     refine int         Level of grid refinement to use. refine = 0 (no grid
                                                                  refinement)
     scheme string      Choose between 'staggered' or 'collocated'
     vtk     int        How often to write vtk files

Generate input data for stencil-based computation.

"""
if length(ARGS) != 5
        println(doc)
        exit(-1)
end

builddir=string(ARGS[1], "/")
outputdir=string(ARGS[2], "/")
refine=parse(Int32, ARGS[3])
scheme=ARGS[4]
vtk_stride=parse(Int32, ARGS[5])

if scheme == "collocated"
        scheme_num = 0
elseif scheme == "staggered"
        scheme_num = 1
else
        error("Unknown scheme selected")
end

import sbp

using Printf
using LinearAlgebra
using SparseArrays

num_moment=4
num_smoothness=4
ref=2^refine
mx=64*ref
my=32*ref
nx=mx + 1
ny=my + 1
r1_s=0.5
r2_s=0.5
r1_r=0.4
r2_r=0.4
dt=2.0/mx
num_steps=250*ref
scheme=0
println(dt)

if scheme_num == 1
        println("Using staggered discretization")
        ac = sbp.StaggeredAcoustic
else
        println("Using collocated discretization")
        ac = sbp.CollocatedAcoustic
end

# Number of grid points and grids
if scheme_num == 1
        nxp, nyp = ac.grid_points("p", nx, ny)
        nx1, ny1 = ac.grid_points("v1", nx, ny)
        nx2, ny2 = ac.grid_points("v2", nx, ny)

        r1_p, r2_p = ac.grids("p", nx, ny, pad=false)
        r1_1, r2_1 = ac.grids("v1", nx, ny, pad=false)
        r1_2, r2_2 = ac.grids("v2", nx, ny, pad=false)
        r1, r2 = ac.grids("node", nx, ny, pad=false)

        pr1_p, pr2_p = ac.grids("p", nx, ny, pad=true)
        pr1_1, pr2_1 = ac.grids("v1", nx, ny, pad=true)
        pr1_2, pr2_2 = ac.grids("v2", nx, ny, pad=true)
        pr1, pr2 = ac.grids("node", nx, ny, pad=true)
else
        nxp = nx
        nyp = ny
        nx1 = nx
        ny1 = ny
        nx2 = nx
        ny2 = ny
        r1_p, r2_p = ac.grid(nx, ny)
        r1_1, r2_1 = ac.grid(nx, ny)
        r1_2, r2_2 = ac.grid(nx, ny)
        r1, r2 = ac.grid(nx, ny)
end

function mapping(r1, r2)
        x = 10.0 .* r1
        sigma = 10/1.05
        y = (5 .+ 1 .* exp.(-sigma^2*(r1 .- 0.5).^2)) .* r2
        return x, y
end

function initial_conditions(xp, yp, x1, y1, x2, y2)
        p = zeros(size(xp))
        v1 = zeros(size(x1))
        v2 = zeros(size(x2))
        return p, v1, v2
end

function zero_pad_x(a, qx, qy)
        b = zeros(qx* qy)

        for i=1:(qx-1)
        for j=1:qy
                b[j + (i - 1) * qy] = a[j + (i - 1) * qy]
        end
        end
        return b
end

function zero_pad_y(a, qx, qy)
        b = zeros(qx* qy)

        for i=1:qx
        for j=1:(qy-1)
                b[j + (i - 1) * qy] = a[j + (i - 1) * (qy - 1)]
        end
        end
        return b
end

x, y = mapping(r1, r2)
xp, yp = mapping(r1_p, r2_p)
x1, y1 = mapping(r1_1, r2_1)
x2, y2 = mapping(r1_2, r2_2)

if scheme_num == 1
        pxp, pyp = mapping(pr1_p, pr2_p)
        px1, py1 = mapping(pr1_1, pr2_1)
        px2, py2 = mapping(pr1_2, pr2_2)
        p, v1, v2 = initial_conditions(pxp, pyp, px1, py1, px2, py2)
        x1 = zero_pad_x(x1, nx + 1, ny + 1)
        y1 = zero_pad_x(y1, nx + 1, ny + 1)
        x2 = zero_pad_y(x2, nx + 1, ny + 1)
        y2 = zero_pad_y(y2, nx + 1, ny + 1)
        
else
        p, v1, v2 = initial_conditions(xp, yp, x1, y1, x2, y2)
end

xs, ys = mapping(r1_s, r2_s)
xr, yr = mapping(r1_r, r2_r)
@printf("Source location: (%f, %f) \n", xs, ys)
@printf("Receiver location: (%f, %f) \n", xr, yr)


sizes = [0, size(p,1), size(v1,1), size(v2,1)]
ranges = accumulate(+, sizes)

@time begin
println("Discretizing in space...")
@printf("refine = %d mx = %d my = %d dt = %g \n", refine, mx, my, dt)
if scheme_num == 1
        ops = ac.init_operators(nx, ny, sbp.OP2019.build_operators, x, y)
else
        st = sbp.Strand1994
        desc_D, desc_H = sbp.Strand1994.load_description(order=4)
        ops = ac.init_operators(nx, ny, st, desc_D, desc_H, x, y)
end
end


"""
Export to C++ code

"""
h1 = 1 / mx
h2 = 1 / my
cfg = Dict("nx" => nxp,
           "ny" => nyp,
           "dt" => dt,
           "hi1" => 1/h1,
           "hi2" => 1/h2,
           "nvtk" => vtk_stride,
           "nt" => num_steps,
           "outputdir" => outputdir,
           "scheme" => scheme_num
                   )

if scheme_num == 0
        r1, h1 = sbp.Grid.grid_xp(nx)
        r2, h2 = sbp.Grid.grid_xp(ny)
        else
        r1, h1 = sbp.Grid.grid_xm(nx + 1)
        r2, h2 = sbp.Grid.grid_xm(ny + 1)
end
d = Array(sbp.Source.source_discretize_2d(r1_s, r2_s, num_moment,
                                    num_smoothness, r1, r2, h1,
                                    h2))
file = x -> string(builddir, x)
        
@time begin
println("Writing discretization to ", builddir)
sbp.lmv.write_config(file("config.txt"), cfg, skip_types=true)

arr = x -> convert(Array{Float32}, Array(diag(x)))
conv = x -> convert(Array{Float32}, Array(x))

if scheme_num == 0
        sbp.lmv.write_vector(file("J.bin"), arr(ops.J))
        sbp.lmv.write_vector(file("g11.bin"), arr(ops.G.b11))
        sbp.lmv.write_vector(file("g12.bin"), arr(ops.G.b12))
        sbp.lmv.write_vector(file("g22.bin"), arr(ops.G.b22))
else
        Jp = arr(ops.Jp)
        J1 = conv(zero_pad_x(arr(ops.J1), nx+1, ny+1))
        J2 = conv(zero_pad_y(arr(ops.J2), nx+1, ny+1))
        g11 = conv(zero_pad_x(arr(ops.G1.b11), nx+1, ny+1))
        g12 = arr(ops.Gp.b12)
        g22 = conv(zero_pad_y(arr(ops.G2.b22), nx+1, ny+1))
        
        sbp.lmv.write_vector(file("Jp.bin"), Jp)
        sbp.lmv.write_vector(file("J1.bin"), J1)
        sbp.lmv.write_vector(file("J2.bin"), J2)
        sbp.lmv.write_vector(file("g1_11.bin"), g11)
        sbp.lmv.write_vector(file("gp_12.bin"), g12)
        sbp.lmv.write_vector(file("g2_22.bin"), g22)
end
        sbp.lmv.write_vector(file("xp.bin"), conv(xp))
        sbp.lmv.write_vector(file("yp.bin"), conv(yp))
        sbp.lmv.write_vector(file("d.bin"), conv(d))
end

sbp.lmv.write_vector(file("v1.bin"), conv(v1))
sbp.lmv.write_vector(file("v2.bin"), conv(v2))
sbp.lmv.write_vector(file("p.bin"), conv(p))
