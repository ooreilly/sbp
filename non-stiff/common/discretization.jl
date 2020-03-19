module Discretization

using LinearAlgebra
using SparseArrays
using DataStructures

import sbp

struct DiscretizationType

        sizes::OrderedDict
        boundary_points::Array
        boundary_conditions::Array
        normal::Array
        boundaries::Array
        mat_boundaries::Array

        function DiscretizationType(names, sizes, boundary_points,
                                           boundary_conditions, normal)

        num_bnd = length(boundary_points)
        @assert num_bnd == length(boundary_conditions)

        vars = variables(names, sizes)

        for indices in boundary_points
        end

        end
end

struct SettingsType
        nx
        dt
        order
        c
        diss
        nu
        nt

        function SettingsType(nx, nu, order, c, T, diss)
                dt = nu / (nx - 1)
                nt = Int64(floor(T/dt))
                return new(nx, dt, order, c, diss, nu, nt)
        end
end

struct OperatorsType

        Dx::AbstractArray
        x::Array
        h::Float64
        Hx::AbstractArray
        Hxi::AbstractArray
        nx::Int64

        function OperatorsType(nx, order)
                Dx, x, h, Hx, Hxi = sbp_operators(nx, order)
                new(Dx, x, h, Hx, Hxi, nx)
        end

        function OperatorsType(nx, order, grid)
                Dx, x, h, Hx, Hxi, mx = staggered_sbp_operators(nx, order, grid)
                return new(Dx, x, h, Hx, Hxi, mx)
        end

end

function normal(nx)
        nl = spzeros(nx)
        nl[1] = -1
        nl[nx] = 1
        return nl
end

function variables(names, sizes)
        d = OrderedDict()

        for (name, size) in zip(names, sizes)
                d[name] = size
        end

        return d
end

function group(vars)
        d = OrderedDict()

        for var in vars
                for (k, v) in var
                        d[k] = v
                end
        end

        return d
end

function boundary(vars::OrderedDict, boundary_points)
        bnd = OrderedDict()
        for (k, v) in vars
                bnd[k] = boundary_points
        end
        return block_diag(restrict_all(vars, bnd))
end

function boundary_each(vars::OrderedDict, boundary_points)
        bnd = OrderedDict()
        i = 1
        for (k, v) in vars
                bnd[k] = boundary_points[i]
                i += 1
        end
        return block_diag(restrict_all(vars, bnd))
end

function restrict_all(vars::OrderedDict, indices::OrderedDict)
        d = OrderedDict()
        
        for (k, v) in vars
                d[k] = restrict(v, indices[k])
        end
        
        return d
end

function restrict(vsize, indices)
        res = spzeros(length(indices), vsize)
        for (j, idx) in enumerate(indices)
                res[j, idx] = 1
        end
        return res
end

function block_diag(vars::OrderedDict)

        rows = zeros(Int64, length(vars))
        cols = zeros(Int64, length(vars))

        i = 1
        for (k, v) in vars
                rows[i] = size(v, 1)
                cols[i] = size(v, 2)
                i += 1
        end

        B = sbp.Sparse.block_matrix(rows, cols)
        
        i = 1
        for (k, v) in vars
                B = sbp.Sparse.block_matrix_insert(B,rows,cols,i,i,v)
                i += 1
        end

        return B
end

function projection_matrix(A, R1, dR1)
        bnd = sbp.BoundaryOperator
        Xp, Lp = bnd.eigen_positive(A)
        Xm, Lm = bnd.eigen_negative(A)
        return projection_matrix_eig(Xp, Xm, R1, dR1)
end

function projection_matrix_eig(Xp, Xm, R1, dR1)
        bnd = sbp.BoundaryOperator
        L = bnd.boundary_operator(Xp, Xm, R1)
        #dLp = bnd.virtual_complementary_boundary_operator(Xp, Xm, dR1)
        dLp = bnd.boundary_operator(-Xp, Xm, dR1)
        P = bnd.virtual_solution_matrix(dLp, L)
        return P
end

function check_energy_stability(H, M)
        Q = H * M + (H * M)'
        return real(eigen(Q).values)
end

function check_stability(M)
        return real(eigen(M).values)
end

function sbp_operators(nx, order)
        # SBP operators
        st=sbp.Strand1994
        desc_D, desc_H = st.load_description(order=order)
        x, h, Dx, Hx = st.build_operators(desc_D, desc_H, nx)
        Hxi = spdiagm(0 => 1 ./ diag(Hx))
        return Dx, x, h, Hx, Hxi
end

function staggered_sbp_operators(nx, order, grid)
        # staggered SBP operators
        st=sbp.OP2019
        desc = st.load_description()
        xp, xm, h, Dp, Dm, Hp, Hm, Pp, Pm = st.build_operators(desc, nx, true)
        if grid == "p"
                Hxi = spdiagm(0 => 1 ./ diag(Hp))
                return Dp, xp, h, Hp, Hxi, nx
        elseif grid == "m"
                Hxi = spdiagm(0 => 1 ./ diag(Hm))
                return Dm, xm, h, Hm, Hxi, (nx + 1)
        else
                error("Unknown grid type")
        end
end

function sat(Lift, A, B, ops, P)
        sat_term_data = kron(Lift, ops.Hxi) *  B' * (A * P)
        sat_term = sat_term_data * B
        return sat_term, sat_term_data
end

function discretization(Ax, ops::OperatorsType)
        M = kron(Ax, -ops.Dx)
        return M
end

function refine(nx)
        return 2 * (nx - 1) + 1
end

function time_step(sys, q, dt, nt, functor=nothing)
        dq = 0 * q
        t = 0.0
        
        for step in 1:nt
                q = sbp.LSRK4.lsrk4(sys, q, dq, t, dt) 
                t = step * dt

                if functor != nothing
                        functor(t, step, q)
                end
        end

        return q, t
end

function time_step_staggered(f, g, u, v, dt, nt::Int64, functor=nothing)
        t = 0.0
        
        for step in 1:nt
                u, v = sbp.StaggeredRK4.staggered_rk4(f, g, u, v, t, dt) 
                t = step * dt
                if functor != nothing
                        functor(t, step, u, v)
                end
        end

        return u, v, t
end

function time_step_staggered_linear(Ap, Av, p, v, dt, nt::Int64, functor=nothing)
        t = 0.0
        
        for step in 1:nt
                p, v = sbp.StaggeredRK4.staggered_rk4_linear(Ap, Av, p, v, dt) 
                t = step * dt
                if functor != nothing
                        functor(t, step, p, v)
                end
        end

        return p, v, t
end

function functional(u, norm)
        return sum(norm * u)
end

function dissipation(A, P)
        diss = P' * A * P
        diss_eig = maximum(abs.(eigen(diss).values))
        return diss_eig
end

function spectral_radius(M, h)
        rho = h * maximum(abs.(eigen(M).values))
        return rho
end
 
function convergence_study(nx0, nt0, dt0, nrefine, errors, rates, init_fcn, disc_fcn, err_fcn, info_fcn=nothing)
        err = []
        nx = nx0
        nt = nt0
        dt = dt0
        grid_points = []
        
        for i in 1:nrefine
                f, ops = disc_fcn(nx)
                u = init_fcn(ops, nx)
                u, t = time_step(f, u, dt, nt)
                err_fcn(errors, u, ops, t)

                if info_fcn != nothing
                        info_fcn(i, nx, nt, dt, ops, errors)
                end
                append!(grid_points, nx)

                nx = refine(nx)
                nt = refine(nt)
                dt = 0.5 * dt
        end
    
        for (k, v) in errors
            rates[k] = sbp.MMS.log2_convergence_rates(v)
        end

        return grid_points
end

function convergence_study_staggered(nx0, nt0, dt0, nrefine, errors, rates, init_fcn, disc_fcn, err_fcn, info_fcn=nothing)
        err = []
        nx = nx0
        nt = nt0
        dt = dt0
        grid_points = []
        
        for i in 1:nrefine
                f, g, ops = disc_fcn(nx)
                u, v = init_fcn(ops, nx)
                u, v, t = time_step_staggered(f, g, u, dt, nt)
                err_fcn(errors, u, v, ops, t)

                if info_fcn != nothing
                        info_fcn(i, nx, nt, dt, ops, errors)
                end
                append!(grid_points, nx)

                nx = refine(nx)
                nt = refine(nt)
                dt = 0.5 * dt
        end
    
        for (k, v) in errors
            rates[k] = sbp.MMS.log2_convergence_rates(v)
        end

        return grid_points
end

end
