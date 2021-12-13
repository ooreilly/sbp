module Bound
using LinearAlgebra


# Estimate constant in Proposition 6 using Frobenius norm.
function estimate_constants(M0, ops)
        c0 = ops.h * LinearAlgebra.opnorm(M0)
        c1 = ops.h * ops.Hxi[1,1]
        return c0, c1
end

end
