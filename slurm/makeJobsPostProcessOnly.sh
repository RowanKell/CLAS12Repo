#!/bin/bash

workdir="/work/clas12/users/gmat/CLAS12Analysis"
#filedir=/work/clas12/users/gmat/CLAS12Analysis/data/fall2018-torus-1-nobg
filedir=/work/clas12/users/gmat/CLAS12Analysis/data/fall2018-torus-1-v1-nSidis
outputdir=${filedir}
rootname="may24_"
processdir="/work/clas12/users/gmat/CLAS12Analysis/macros/dihadron_process/"
processcodename="pipi0_postprocess_only.C"

runJobs="${workdir}/slurm/runJobsPostProcessOnly.sh"
touch $runJobs
chmod +x $runJobs
echo " " > $runJobs

i=0

for rootfile in "$filedir"/may24_*
do
    file="${workdir}/slurm/shells/${rootname}${i}.sh"
    touch $file
    echo "#!/bin/tcsh" > $file
    echo "#SBATCH --account=clas12" >> $file
    echo "#SBATCH --partition=production" >> $file
    echo "#SBATCH --mem-per-cpu=4000" >> $file
    echo "#SBATCH --job-name=${rootname}${i}" >> $file
    echo "#SBATCH --cpus-per-task=4" >> $file
    echo "#SBATCH --time=24:00:00" >> $file
    echo "#SBATCH --chdir=${workdir}" >> $file
    echo "#SBATCH --output=${workdir}/slurm/output/%x-%j-%N.out" >> $file
    echo "#SBATCH --error=${workdir}/slurm/output/%x-%j-%N.err" >> $file
    echo "echo ${workdir}" >> $file
    echo "source /group/clas12/packages/setup.csh" >> $file
    echo "module load clas12/pro" >> $file
    echo "set CLAS12ROOT=/w/hallb-scshelf2102/clas12/users/gmat/packages/clas12root" >> $file
    echo "set CCDB_HOME=${CLAS12ROOT}/ccdb" >> $file
    echo "source ${CCDB_HOME}/environment.csh" >> $file
    echo "cd ${processdir}" >> $file    
    echo "clas12root ${processcodename}\\(\\\"${rootfile}\\\"\\)" >> $file   
    echo "sbatch shells/${rootname}${i}.sh" >> $runJobs
    i=$((i+1))
done
