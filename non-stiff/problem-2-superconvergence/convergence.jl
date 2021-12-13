import sbp
using MAT
using LinearAlgebra
include("../common/wave_equation.jl")
include("../common/bound.jl")
include("../common/discretization.jl")
disc = Discretization
pde=WaveEquation

datafile=ARGS[1]
dR = parse(Float64, ARGS[2]) * ones(1,1)
nsamples = 4
nrefinements = 6
nx = 20
T = 1.2
dt = 0.25 / (nx - 1)
nt = T / dt
order = 4

mms = pde.manufactured_solution

function data()
        return Dict("p" => [], "v" => [] , 
                    "Fp" => [], "Fv" => [])
end

function discretize(nx, dRi)
    ops = disc.OperatorsType(nx, order)
    M = disc.discretization(pde.A(), ops)
    Rl = -ones(1,1)
    Rr = -ones(1,1)
    Pl = disc.projection_matrix(pde.A(pde.LEFT), Rl, dRi)
    Pr = disc.projection_matrix(pde.A(pde.RIGHT), Rr, dRi)
    Sl, Slg = pde.sat(pde.A(pde.LEFT), pde.LEFT, ops, Pl)
    Sr, Srg = pde.sat(pde.A(pde.RIGHT), pde.RIGHT, ops, Pr)
    gl = t -> pack(mms(ops.x[1], t))
    gr = t -> pack(mms(ops.x[end], t))
    f = (u,t) -> (M + Sl + Sr) * u - Slg * gl(t) - Srg * gr(t)
    return f, ops
end

function pack(mms_)
   p = mms_["p"] 
   v = mms_["v"] 
   return [p; v]
end
        

function initialize(ops, nx)
   return pack(mms(ops.x, 0))
end

function error(d, u, ops, t)
    mms = pde.manufactured_solution(ops.x,  t)
    F_mms = pde.functionals(t)

    soln = Dict(sbp.MMS.extract_vars(u, pde.sizes(ops.nx)))
    norms = Dict("p" => ops.Hx, "v" => ops.Hx)
    err = sbp.MMS.l2_errors(soln, mms, norms)
    F_err = sbp.MMS.functional_errors(soln, F_mms, norms)
    append!(d["p"], err["p"])
    append!(d["v"], err["v"])
    append!(d["Fp"], F_err["p"])
    append!(d["Fv"], F_err["v"])
end

function info(i, nx, nt, dt, ops, err)
        println("Grid = $i, nx = $nx, p = ", err["p"][end],
                " v = ", err["v"][end], 
                ", rate p = ", sbp.MMS.log2_convergence_rates(err["p"])[end], 
                ", rate v = ", sbp.MMS.log2_convergence_rates(err["v"])[end],
                ", rate Fp = ", sbp.MMS.log2_convergence_rates(err["Fp"])[end], 
                ", rate Fv = ", sbp.MMS.log2_convergence_rates(err["Fv"])[end]
                )
end

errors = data()
rates = data()

x = disc.convergence_study(nx, nt, dt, nrefinements, errors, rates,
                           initialize, x -> discretize(x, dR), error,
                           info)
println("")
matwrite(datafile, 
         Dict("dR" => dR, 
              "x" => x,
              "p" => errors["p"], 
              "v" => errors["v"], 
              "Fp" => errors["Fp"], 
              "Fv" => errors["Fv"], 
              "q_p" => rates["p"],
              "q_v" => rates["v"],
              "q_Fp" => rates["Fp"],
              "q_Fv" => rates["Fv"]
               ))

