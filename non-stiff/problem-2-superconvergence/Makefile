JL=julia
PY=python3     
show_fig=0
dR_dissipative=0
dR_conservative=1
dissipative_data="data/dissipative.mat"
conservative_data="data/conservative.mat"

all: init dissipative dissipative-table conservative conservative-table

define crop
	pdfcrop $(1) $(1)
endef

init:
	mkdir -p figures data

dissipative:
	$(JL) convergence.jl $(dissipative_data) $(dR_dissipative)

dissipative-table:
	$(JL) convergence_table.jl $(dissipative_data)

conservative:
	$(JL) convergence.jl $(conservative_data) $(dR_conservative)

conservative-table:
	$(JL) convergence_table.jl $(conservative_data)

clean:
	rm data/*
