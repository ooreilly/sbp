for ((i = 0; i < 7; i++))
do
zip results.zip timing/collocated/output_default_$i/log.txt
zip results.zip timing/collocated/output_default_$i/u.csv
zip results.zip timing/staggered/output_modified_$i/log.txt
zip results.zip timing/staggered/output_modified_$i/u.csv
zip results.zip timing_stencil/staggered/output_modified_$i/log.txt
done
