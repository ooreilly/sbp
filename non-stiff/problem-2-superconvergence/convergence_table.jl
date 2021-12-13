using MAT
import sbp 
datafile=ARGS[1]
 
vars = matread(datafile)

nx = vars["x"]
p = vars["p"]
v = vars["v"]
qp = vars["q_p"]
qv = vars["q_v"]
Fp = vars["Fp"]
Fv = vars["Fv"]
qFp = vars["q_Fp"]
qFv = vars["q_Fv"]

header = [ "Nx", "p", "v",   
           "rate p", "rate v",
          "Fp", "rate Fp", 
          "Fv", "rate Fv"]
variables = [
             nx, 
             round.(log10.(p), sigdigits=4), 
             round.(qp, sigdigits=4),
             round.(log10.(v), sigdigits=4), 
             round.(qv, sigdigits=4),
             round.(log10.(Fp), sigdigits=4), 
             round.(qFp, sigdigits=4),
             round.(log10.(Fv), sigdigits=4), 
             round.(qFv, sigdigits=4),
             ]
print(join(sbp.Table.latex(header, variables)))
