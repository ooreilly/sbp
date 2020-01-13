for ((i=6; i < 7; i++))
do
make refine=$i task=timing vtk_stride=0 scheme=staggered metric_tensor=modified
make refine=$i task=timing vtk_stride=0 scheme=staggered metric_tensor=default
make refine=$i task=timing vtk_stride=0 scheme=collocated metric_tensor=default
done
