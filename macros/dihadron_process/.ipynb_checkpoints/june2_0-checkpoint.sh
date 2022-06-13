#!/bin/tcsh
echo /work/clas12/users/rojokell/CLAS12Analysis
source /group/clas12/packages/setup.csh
module load clas12/pro
set CLAS12ROOT=/w/hallb-scshelf2102/clas12/users/rojokell/clas12root
set CCDB_HOME=/group/clas12/packages/clas12root/1.7.8.b/ccdb
#source /site/12gev_phys/devel/Linux_CentOS7.7.1908-x86_64-gcc9.2.0/ccdb/ccdb-1.07.00/environment.csh
#cd /work/clas12/users/rojokell/CLAS12Analysis/macros/dihadron_process/
clas12root rowanTestpipluspiminus.C\(\"/cache/clas12/rg-a/production/recon/fall2018/torus+1/pass1/v1/dst/train/nSidis//nSidis_005422.hipo\",\"/work/clas12/users/rojokell/CLAS12Analysis/data/test/june2_0.root\"\)
