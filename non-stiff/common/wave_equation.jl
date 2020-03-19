module WaveEquation
include("discretization.jl")
using LinearAlgebra
disc=Discretization

names = ["p", "v"]
WAVE_NUMBER = 8 * pi
LEFT = -1
RIGHT = 1



function A(n)
        Aout = zeros(2,2)
        Aout[1,2] = n
        Aout[2,1] = n
        return Aout
end

function A()
        A = zeros(2,2)
        A[1,2] = 1
        A[2,1] = 1
        return A
end

function reflection_matrix(r)
        R1 = zeros(1,1)
        R1[1,1] = r
        return R1
end

function projection(a, dr::Float64)
        dr_ = ones(1,1) * dr
        r = -ones(1,1)
        P = projection_matrix(A(a), r, dr_)
        return P
end

function sat_left(ops::disc.OperatorsType, P::Matrix)
        B = boundary(sizes(ops.nx), [1])
        I2 = Matrix(I, 2, 2)
        sat_term_data = kron(I2, ops.Hxi) * B' * (A(LEFT) * P)
        sat_term = sat_term_data * B
        return sat_term, sat_term_data
end

function sat_right(ops::disc.OperatorsType, P::Matrix)
        B = boundary(sizes(ops.nx), [ops.nx])
        I2 = Matrix(I, 2, 2)
        sat_term_data = kron(I2, ops.Hxi) * B' * (A(RIGHT) * P)
        sat_term = sat_term_data * B
        return sat_term, sat_term_data
end

function sat(A, side, ops, P)
        if side == LEFT
                B = disc.boundary(sizes(ops.nx), [1])
        else
                B = disc.boundary(sizes(ops.nx), [ops.nx])
        end
        I_ = Matrix(I, 2, 2)
        return disc.sat(I_, A, B, ops, P)
end

function zero_p_left()
        return reflection_matrix(-1)
end

function zero_v_right()
        return reflection_matrix(-1)
end

function non_reflective()
        return zeros(1,1)
end

function leakage(alpha)
        return ones(1,1) * (alpha - 1) / (1 + alpha)
end

function naive_leakage(alpha, ops)
        B = disc.boundary(sizes(ops.nx), [ops.nx])
        P = zeros(2, 2)
        P[2,1] = 1 
        P[2,2] = -alpha
        I_ = Matrix(I, 2, 2)
        s1, s2 = disc.sat(I_, I_, B, ops, P)
        return P, s1, s2 
end

function energy_norm(Hx)
        return kron(Matrix(I, 2, 2), Hx)
end

function manufactured_solution(x, t)
        k = WAVE_NUMBER
        p = t -> cos.(k * t) .* sin.(k * x)
        v = t -> - sin.(k * t) .* cos.(k * x)
        q = Dict("p" => p(t), "v" => v(t))
        return q
end

function functionals(t)
        k = WAVE_NUMBER
        Fp = (1 - cos(k)) / k * cos(k * t)
        Fv = - sin(k) / k * sin(k * t)
        return Dict("p" => Fp, "v" => Fv)
end

function sizes(nx)
        return disc.variables(names, [nx, nx])
end

function restrict_left(nx)
        left = disc.variables(names, [1, 1])
        block_left = disc.block_diag(restrict(vars, left))
end

end
