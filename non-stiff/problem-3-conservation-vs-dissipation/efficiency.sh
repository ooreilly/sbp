#for i in  `seq 0 5`
for i in  `seq 4 5`
do
        make rk4_ed refine=$i | tee logs/rk4_ed_$i.txt
        #make srk4_ec refine=$i | tee logs/srk4_ec_$i.txt
        #make srk4_opt refine=$i | tee logs/srk4_opt_$i.txt
        #make rk4_naive refine=$i | tee logs/rk4_naive_$i.txt
done
