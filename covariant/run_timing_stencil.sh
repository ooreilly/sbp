for ((i=0; i < 7; i++))
do
make refine=$i task=timing_stencil  vtk_stride=0 scheme=staggered init build-stencil run-stencil
make refine=$i task=timing_stencil  vtk_stride=0 scheme=collocated init build-stencil run-stencil
done
