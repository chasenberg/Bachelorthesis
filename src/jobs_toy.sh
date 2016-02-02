#!/bin/sh

#PBS -N %(job_name)s
#PBS -o %(out_file)s
#PBS -e %(err_file)s
#PBS -m n
#PBS -M phidomuell@googlemail.com
#PBS -j oe
#PBS -l walltime=%(walltime)s
#PBS -l nodes=1:ppn=%(num_cpu)s
#PBS -l vmem=%(vmem)s
export SCRATCH=/local/$USER/$PBS_JOBID

echo "Auto-generated PBS script"
echo "B2Dpi FT toy studies starting at `date`, running on node `hostname`"
source /lhcbsoft/tu-dortmund/LHCbSoftwareSetup.sh
time . SetupProject.sh DaVinci v36r1
echo $0
echo $PBS_JOBID > %(log_file)s
echo ${PBS_JOBID} >> %(log_file)s	
cd %(cwd)s

for i in `seq %(seeds)s`;
do
  time /./home/chasenberg/Repository/bachelor-template/build/bin/dootoycp_spline -c %(jobs_dir)s/dootoycp-config.txt --toyfac.random_seed=$i >> %(log_file)s
done

echo finished job number %(job_number)s
