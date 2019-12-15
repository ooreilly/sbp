refine=0
task=default
scheme=staggered
metric_tensor=modified
builddir=$(task)/$(scheme)/build_$(metric_tensor)_$(refine)
outputdir=$(task)/$(scheme)/output_$(metric_tensor)_$(refine)
run_w_julia=0
vtk_stride=10

all: init build run

init:
	mkdir -p $(builddir)
	mkdir -p $(outputdir)

build:
	julia simulation.jl $(builddir) $(outputdir) $(refine) $(scheme) \
		$(metric_tensor) $(run_w_julia) $(vtk_stride)

run:
	build/src/solver $(builddir)

.PHONY: build
