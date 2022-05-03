Instructions for installation of CLAS12Analysis software

1. Ensure both ROOTSYS and CLAS12ROOT are environment variables (see https://github.com/JeffersonLab/clas12root for details on installing CLAS12ROOT). For ROOT, we need 6.24.06, which should automatically install with 'module load clas12/pro' (No need to source /apps/root/6.22.06/bin/thisroot.csh)

2. 'cd' into the build directory

3. Run ./autogen.sh --prefix=/path/to/CLAS12Analysis/

4. Run make

5. Run make install

6. 'cd ..' and ensure both a lib/ and include/ directory has been created

7. In your home directory, edit the .cshrc, add the following line
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/path/to/CLAS12Analysis/lib

8. Reopen your terminal or source ~/.cshrc to have the LD_LIBRARY_PATH appended to
