doc = 
"""
usage: 
simulation.jl <outputdir> <builddir> <refine> <scheme> <metric_tensor> <run> <vtk>

     builddir path      Directory to write discretization.
     outputdir path     Directory to write output data to.
     refine int         Level of grid refinement to use. refine = 0 (no grid
                                                                  refinement)
     scheme string      Choose between 'staggered' or 'collocated'
     metric_tensor string  Set to 'modified' to use the modified metric tensor
                           discretization.      
     run     int        1 (Run simulation in Julia) 0 (Otherwise) 
     vtk     int        How often to write vtk files

This example code solves the acoustic wave equation in covariant form on a
curvilinear staggered grid

The governing equations are discretized in space and the zero pressure boundary
condition is applied on all boundaries. The resulting discretization is stored
in a sparse matrix `A`. The solution is advanced in time by solving the linear
system

dq/dt = A*q

using a low-storage Runge-Kutta time integration scheme.
"""

"""
Parse command line arguments
"""

if length(ARGS) != 7
        println(doc)
        exit(-1)
end

builddir=string(ARGS[1], "/")
outputdir=string(ARGS[2], "/")
refine=parse(Int32, ARGS[3])
scheme=ARGS[4]
metric_tensor=ARGS[5]
run=parse(Int32, ARGS[6])
vtk_stride=parse(Int32, ARGS[7])

"""
User settings

nx, ny: Number of grid points in each direction (nodes)
dt: Time step size
num_steps: Number of time steps to take
info_stride: How often to display information
outputdir: Save vtk snapshots to outputdir/{p, v1, v2}_%d.vtk
r1_s, r2_s: Source location in parameter space
r1_r, r2_r: Receiver location in parameter space
num_moment: Number of moment conditions used in source discretization
num_smoothness Number of smoothness conditions used in source discretization
run: Run simulation in Julia
export: Export the discretization for running with the C++ code
"""

ref=2^refine
mx=64*ref
my=32*ref
nx=mx + 1
ny=my + 1
dt=2/mx
num_steps=250*ref
info_stride=100*ref
vtk_stride=vtk_stride*ref
num_moment=4
num_smoothness=4
r1_s=0.5
r2_s=0.5
r1_r=0.4
r2_r=0.4
export_all=true

"""
Define the mapping function x = X(r^1, r^2), y = Y(r^1, r^2) that describes the
geometry of the computational domain.
"""
function mapping(r1, r2)
        x = 10.0 .* r1
        sigma = 10/1.05
        y = (5 .+ 1 .* exp.(-sigma^2*(r1 .- 0.5).^2)) .* r2
        return x, y
end

"""
Define the initial conditions for p, v1, v2. Make sure that p is of the same
size as the grid coordinate vectors (xp, yp), and similarly for v1, v2.
"""                                 
function initial_conditions(xp, yp, x1, y1, x2, y2)
        p = zeros(size(r1_p))
        v1 = zeros(size(r1_1))
        v2 = zeros(size(r1_2))
        return p, v1, v2
end

"""
Define the source time function that introduces a time dependent forcing in the
pressure equation via a singular source term.
"""
function source_time_function(t)
        return exp.( -5e2 .* (t  .- 2.0).^2)
end

###############################################################################

using Printf
using SparseArrays
import sbp
using sbp.Grid: grid_xp, grid_2d_x, grid_2d_y
using sbp.VTK: vtk_write
using sbp.LSRK4:lsrk4

if scheme == "staggered"
        println("Using staggered discretization")
        ac = sbp.StaggeredAcoustic
else
        println("Using collocated discretization")
        ac = sbp.CollocatedAcoustic
end

# Number of grid points and grids
if scheme == "staggered"
        nxp, nyp = ac.grid_points("p", nx, ny)
        nx1, ny1 = ac.grid_points("v1", nx, ny)
        nx2, ny2 = ac.grid_points("v2", nx, ny)
        r1_p, r2_p = ac.grids("p", nx, ny)
        r1_1, r2_1 = ac.grids("v1", nx, ny)
        r1_2, r2_2 = ac.grids("v2", nx, ny)
        r1, r2 = ac.grids("node", nx, ny)
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

x, y = mapping(r1, r2)
xp, yp = mapping(r1_p, r2_p)
x1, y1 = mapping(r1_1, r2_1)
x2, y2 = mapping(r1_2, r2_2)

# Pack initial conditions into `q = [p, v1, v2]`
p, v1, v2 = initial_conditions(xp, yp, x1, y1, x2, y2)
q = vcat(p, v1, v2)
sizes = [0, size(p,1), size(v1,1), size(v2,1)]
ranges = accumulate(+, sizes)

"""

Discretize in space using either the staggered or collocated scheme

"""

@time begin
println("Discretizing in space...")
@printf("refine = %d mx = %d my = %d dt = %g \n", refine, mx, my, dt)
if scheme == "staggered"
        ops = ac.init_operators(nx, ny, sbp.OP2019.build_operators, x, y)
else
        st = sbp.Strand1994
        desc_D, desc_H = sbp.Strand1994.load_description(order=4)
        ops = ac.init_operators(nx, ny, st, desc_D, desc_H, x, y)
end
if metric_tensor == "modified" && scheme == "staggered"
        println("using modified metric tensor discretization")
        Gp = ac.modified_contravariant_metric_tensor(ops)
else
        println("using positivity preserving metric tensor discretization")
        Gp = ac.contravariant_metric_tensor(ops)
end

Ap = ac.divergence(ops)
Av = ac.gradient(ops, Gp)
S = ac.pressure_sat(ops, Gp)
A = -ac.spatial_discretization(Ap, Av - S)
d = ac.pressure_source(ops, r1_s, r2_s, num_moment, num_smoothness)
r = ac.pressure_receiver(ops, r1_r, r2_r, num_moment, num_smoothness)
xs, ys = mapping(r1_s, r2_s)
xr, yr = mapping(r1_r, r2_r)
@printf("Source location: (%f, %f) \n", xs, ys)
@printf("Receiver location: (%f, %f) \n", xr, yr)

end

A = dropzeros(A)

println("Number of non-zeros in A: ", length(A.nzval))

"""
Export to C++ code

"""

if export_all == true
        cfg = Dict("dt" => dt, 
                   "nxp" => nxp,
                   "nyp" => nyp,
                   "ninfo" => 20*ref,
                   "ncsv" => ref,
                   "nvtk" => vtk_stride,
                   "ndump" => num_steps,
                   "nrecv" => 1,
                   "nt" => num_steps,
                   "outputdir" => outputdir, 
                   )
        
        println("Writing discretization to ", builddir)
        sbp.lmv.write_config(string(builddir,  "config.txt"), cfg,
                                                               skip_types=true)
        sbp.lmv.write_vector(string(builddir, "q.bin"), q)
        Aout = convert(SparseMatrixCSC{Float64, Int32}, A)
        sbp.lmv.write_csr_matrix(string(builddir, "A.bin"), Aout)
        r = Array(r)[:]
        sbp.lmv.write_vector(string(builddir, "r.bin"), r)
        d = Array(d)[:]
        sbp.lmv.write_vector(string(builddir, "d.bin"), d)
        sbp.lmv.write_vector(string(builddir, "xp.bin"), xp)
        sbp.lmv.write_vector(string(builddir, "yp.bin"), yp)

end

if run == false
        exit()
end

"""
Advance solution in time

"""

scheme = (u, t) -> A * u + d * source_time_function(t)

dq = zeros(size(q))
curr_vtk = 0
elapsed = 0
@printf("step \t time \t time per step \t elapsed time \n")
for step in 1:num_steps 

        tk = (step - 1) * dt

        if step % vtk_stride == 0
                outfile=field -> @sprintf("%s/%s_%d.vtk", outputdir, field,
                                          curr_vtk)
                p = q[1+ranges[1]:ranges[2]]
                v1 = q[1+ranges[2]:ranges[3]]
                v2 = q[1+ranges[3]:ranges[4]]
                vtk_elapsed = @elapsed begin
                        vtk_write(outfile("p"), xp, yp, p, nxp, nyp)
                        vtk_write(outfile("v1"), x1, y1, v1, nx1, ny1)
                        vtk_write(outfile("v2"), x2, y2, v2, nx2, ny2)
                end
                @printf("Wrote: %s in %.3f s \n", outfile("{p, v1, v2}"),
                        vtk_elapsed)
                global curr_vtk += 1
        end

        global elapsed += @elapsed global q = lsrk4(scheme, q, dq, tk, dt)
        if step % info_stride == 0
                @printf("%05d \t %.3f \t %.3f s \t %.3f s \n", step, tk, elapsed
                        / step, elapsed)
        end

end
