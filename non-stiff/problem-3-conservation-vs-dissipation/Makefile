exe=build/acoustic

refine=0
nx=102
ny=202
Lx=1.0
Ly=2.0
T=0.5
vtk_stride=-100
stride=1
nu_vtk=0.4

# Source settings
t0=0.25
fp=6.4
sx=0.45
sy=1.0

# CFL numbers
nu_naive=0.068
nu_rk4_opt = 0.412
nu_srk4_opt = 1.584
nu_srk4_opt = 1.725
nu_srk4_ec3 = 1.725

nu_rk4_ec3 = 1.0



all: init lf srk4 rk4

init:
	mkdir -p output
	mkdir -p logs
	mkdir -p figures

rk4_ed:
	$(exe) output/rk4_ed_$(refine) $(refine) 2 3 $(nx) $(ny) $(Lx) $(Ly) $(nu_rk4_opt) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

rk4_naive:
	$(exe) output/rk4_naive_$(refine) $(refine) 2 4 $(nx) $(ny) $(Lx) $(Ly) $(nu_naive) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

srk4_ec:
	$(exe) output/srk4_ec_$(refine) $(refine) 1 2 $(nx) $(ny) $(Lx) $(Ly) $(nu_srk4_opt) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

srk4_opt:
	$(exe) output/srk4_opt_$(refine) $(refine) 1 5 $(nx) $(ny) $(Lx) $(Ly) $(nu_srk4_opt) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

srk4_ec3:
	$(exe) output/srk4_ec3_$(refine) $(refine) 1 6 $(nx) $(ny) $(Lx) $(Ly) $(nu_srk4_ec3) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

rk4_ec3:
	$(exe) output/rk4_ec3_$(refine) $(refine) 2 6 $(nx) $(ny) $(Lx) $(Ly) $(nu_rk4_ec3) $(vtk_stride) $(stride) $(T) $(fp) $(t0) $(sx) $(sy)


# Save vtk files for showing interface inaccuracy problem
vtk_ec:
	$(exe) output/vtk_ec_$(refine) $(refine) 2 2 $(nx) $(ny) $(Lx) $(Ly) $(nu_vtk) 10 $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

vtk_opt:
	$(exe) output/vtk_opt_$(refine) $(refine) 2 5 $(nx) $(ny) $(Lx) $(Ly) $(nu_vtk) 10 $(stride) $(T) $(fp) $(t0) $(sx) $(sy)

vtk_ed:
	$(exe) output/vtk_ed_$(refine) $(refine) 2 3 $(nx) $(ny) $(Lx) $(Ly) $(nu_vtk) 10 $(stride) $(T) $(fp) $(t0) $(sx) $(sy)


clean:
	rm output/*

