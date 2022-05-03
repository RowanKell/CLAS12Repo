Instructions for installation of CLAS12Analysis software

1. Ensure both ROOTSYS and CLAS12ROOT are environment variables (see https://github.com/JeffersonLab/clas12root for details on installing CLAS12ROOT). For ROOT, we need 6.24.06, which should automatically install with 'module load clas12/pro'

2. 'cd' into the build directory

3. Run ./autogen.sh --prefix=/path/to/CLAS12Analysis/

4. Run make

5. Run make install

6. 'cd ..' and ensure both a lib/ and include/ directory has been created
