Instructions for installation of dihadron analysis software

1. In your home directory, create an empty 'install' folder

2. Ensure both ROOTSYS and CLAS12ROOT are environment variables (see https://github.com/JeffersonLab/clas12root for details on installing CLAS12ROOT). For ROOT, we need 6.24.06, which should automatically install with 'module load clas12/pro'

3. Type the command 'aclocal' and press ENTER (could likely just put this in autogen.sh but I'm lazy ;)

4. Run ./autogen.sh --prefix=/path/to/home/install

5. Run make

6. Run make install

7. Ensure both a lib/ and include/ directory has been created with dihadronana files stored inside
