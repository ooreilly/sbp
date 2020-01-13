i=4
make refine=$i task=visualization scheme=staggered metric_tensor=modified vtk_stride=10
make refine=$i task=visualization scheme=staggered metric_tensor=default  vtk_stride=10
make refine=$i task=visualization scheme=collocated metric_tensor=default vtk_stride=10
