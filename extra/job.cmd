#@ output = job.out
#@ total_tasks = 1
#@ node = 1
#@ wall_clock_limit = 00:15:00
#@ class = micro
#@ job_name = gccg
#@ job_type = parallel
#@ notify_user = miklos.homolya@tum.de
#@ queue
#@ initialdir = /home/hpc/h039v/h039vam/team06/A1
for dir in code_*; do
    cd ${dir}
    ./gccg text tjunc.dat tjunc.run1_
    ./gccg text tjunc.dat tjunc.run2_
    ./gccg text tjunc.dat tjunc.run3_
    ./gccg text drall.dat drall.run1_
    ./gccg text drall.dat drall.run2_
    ./gccg text drall.dat drall.run3_
    ./gccg text pent.dat pent.run1_
    ./gccg text pent.dat pent.run2_
    ./gccg text pent.dat pent.run3_
    ./gccg text cojack.dat cojack.run1_
    ./gccg text cojack.dat cojack.run2_
    ./gccg text cojack.dat cojack.run3_
    cd ..
done
