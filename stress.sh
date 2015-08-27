export OMP_NUM_THREADS=4
for d in `seq 10`;
do 
 /usr/bin/time ./neutral-test --popsize 100000 --numloci 5 --inittraits 6 --innovrate 0.0001 --simlength 25000 --debug 1 --ruletype wfia
 sleep 1
done
