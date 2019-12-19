refine=0
task=default
scheme=collocated
metric_tensor=modified
builddir=$(task)/$(scheme)/build_$(metric_tensor)_$(refine)
outputdir=$(task)/$(scheme)/output_$(metric_tensor)_$(refine)
run_w_julia=0
vtk_stride=10

ifeq ($(scheme), staggered)
stencil_solver=staggered_solver
else
stencil_solver=collocated_solver
endif

all: init build run

init:
	mkdir -p $(builddir)
	mkdir -p $(outputdir)

build:
	julia simulation.jl $(builddir) $(outputdir) $(refine) $(scheme) \
		$(metric_tensor) $(run_w_julia) $(vtk_stride)

run:
	build/src/solver $(builddir)

build-stencil:
	julia stencil.jl $(builddir) $(outputdir) $(refine) $(scheme) \
		$(vtk_stride) 

run-stencil:
	build/src/$(stencil_solver) $(builddir)

.PHONY: build
