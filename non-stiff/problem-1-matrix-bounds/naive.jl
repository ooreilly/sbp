import sbp
using MAT
using LinearAlgebra
include("../common/wave_equation.jl")
include("../common/bound.jl")
include("../common/discretization.jl")
disc = Discretization
pde=WaveEquation

datafile = ARGS[1]
nx = 100
nsamples = 100
alpha = range(0.0, stop=10.0, length=nsamples)
order = 4

Mnorm = zeros(nsamples)
sigma_norm_upper = zeros(nsamples)
sigma_norm_lower = zeros(nsamples)

ops = disc.OperatorsType(nx, order)
h = ops.h

M0 = disc.discretization(pde.A(), ops)

sigma_d = LinearAlgebra.svd(Matrix(ops.Dx) * h)
rd = sigma_d.S[2]
sd = sigma_d.S[1]
rl = ops.Hxi[1,1] * h
sl = rl

# Left penalty term
Pl = disc.projection_matrix(pde.A(pde.LEFT), zeros(1,1), zeros(1,1))
Ml, sat_term_data = pde.sat(pde.A(pde.LEFT), pde.LEFT, ops, Pl)

Sl =  pde.A() * Pl

println("Estimated constants: $rd <= h\\|D\\|_2 <= $sd")
println("Estimated constants: $rl <= h\\|L\\|_2 <= $sl")

i = 1
for alphai in alpha

        # Right penalty term
        Sr, Mr, sat_term_data = pde.naive_leakage(alphai, ops)

        M = M0 + Ml + Mr

        sigma_norm_upper[i] = LinearAlgebra.svd(Matrix(Sr)).S[1] + LinearAlgebra.svd(Matrix(Sl)).S[1] 
        sigma_norm_lower[i] = LinearAlgebra.svd(Matrix(Sr)).S[1] - LinearAlgebra.svd(Matrix(Sl)).S[1]
        
        eig_max = maximum(disc.check_stability(M))
        if (eig_max > 0) 
                println("Warning: Unstable. eig_max = $eig_max")
        end
        Mnorm[i] = h * LinearAlgebra.opnorm(M)
        println("\\|M\\|_2 = ", Mnorm[i], " ", sigma_norm_lower[i], " <= \\|Sigma\\|_2 <= ", sigma_norm_upper[i])
        global i += 1
end
matwrite(datafile, Dict("alpha" => collect(alpha), 
                                   "Mnorm" => Mnorm, 
                                   "rd" => rd,
                                   "sd" => sd,
                                   "rl" => rl,
                                   "sl" => sl,
                                   "sigma_norm_upper" => sigma_norm_upper,
                                   "sigma_norm_lower" => sigma_norm_lower)
         )
