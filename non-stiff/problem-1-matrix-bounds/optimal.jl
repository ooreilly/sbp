import sbp
using MAT
using LinearAlgebra
include("../common/wave_equation.jl")
include("../common/bound.jl")
include("../common/discretization.jl")
disc = Discretization
pde = WaveEquation

datafile = ARGS[1]
nx = 100
nsamples = 100
order = 4
alpha = 10
dR = range(-1.0, stop=1.0, length=nsamples)

vars = pde.sizes(nx)
boundary_points = [1, nx]
boundary_conditions = [pde.non_reflective, () -> pde.leakage(alpha)]
sides = [pde.LEFT, pde.RIGHT]
boundaries = [disc.boundary(vars, [1]), disc.boundary(vars, [nx])]
nl = disc.normal(nx)



ops = disc.OperatorsType(nx, order)
h = ops.h
M0 = disc.discretization(pde.A(), ops)

sigma_d = LinearAlgebra.svd(Matrix(ops.Dx) * h)
rd = sigma_d.S[2]
sd = sigma_d.S[1]
rl = ops.Hxi[1,1] * h
sl = rl
println("estimated constants: $rd <= h\\|d\\|_2 <= $sd")
println("estimated constants: $rl <= h\\|l\\|_2 <= $sl")
i = 1

Mnorm = zeros(nsamples)
sigma_norm_left = zeros(nsamples)
sigma_norm_right = zeros(nsamples)
sigma_norm_upper = zeros(nsamples)
sigma_norm_lower = zeros(nsamples)
for dRi  in dR

        coef = [0, dRi]

        M = M0

        """
         Go through each boundary and construct projection matrices, and sat
         terms
        """
        for (dRi, bnd, bc, bp, side) in zip(coef, boundaries, boundary_conditions,
                                      boundary_points, sides)
                r = bc()
                dr = pde.reflection_matrix(dRi)
                Ai = pde.A(nl[bp])
                P = disc.projection_matrix(Ai, r, dr)
                sat_term, sat_term_data = pde.sat(Ai, side, ops, P)
                M = M + sat_term
                if side == pde.RIGHT
                 sigma_norm_right[i] = LinearAlgebra.opnorm(Ai * P)
                end
                if side == pde.LEFT
                 sigma_norm_left[i] = LinearAlgebra.opnorm(Ai * P)
                end
        end

        sigma_norm_upper[i] = sigma_norm_left[i] + sigma_norm_right[i]
        sigma_norm_lower[i] = sigma_norm_right[i] - sigma_norm_left[i]
        
        # Spectral radius
        eig_max = maximum(disc.check_stability(M))
        if (eig_max > 0) 
                println("Warning: Unstable. eig_max = $eig_max")
        end
        Mnorm[i] = h * LinearAlgebra.opnorm(M)
        println("\\|M\\|_2 = ", Mnorm[i], " ", sigma_norm_lower[i], " <= \\|Sigma\\|_2 <= ", sigma_norm_upper[i])
        global i += 1
end
matwrite(datafile, Dict("dR" => collect(dR), 
                                   "Mnorm" => Mnorm, 
                                   "rd" => rd,
                                   "sd" => sd,
                                   "rl" => rl,
                                   "sl" => sl,
                                   "sigma_norm_upper" => sigma_norm_upper,
                                   "sigma_norm_lower" => sigma_norm_lower)
         )
