JL=julia
PY=python3     
show_fig=0

optimal_data="data/optimal.mat"
optimal_fig="figures/problem1_optimal.pdf"

naive_data="data/naive.mat"
naive_fig="figures/problem1_naive.pdf"

all: init optimal optimal-fig naive naive-fig

define crop
	pdfcrop $(1) $(1)
endef

init:
	mkdir -p data figures

optimal:
	$(JL) optimal.jl $(optimal_data)

optimal-fig:
	$(PY) plot_optimal.py $(optimal_data) $(show_fig) $(optimal_fig)
	$(call crop, $(optimal_fig))

naive:
	$(JL) naive.jl $(naive_data)

naive-fig:
	$(PY) plot_naive.py $(naive_data) $(show_fig) $(naive_fig)
	$(call crop, $(naive_fig))

clean:
	rm figures/* data/*
