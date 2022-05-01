Instructions for installation of dihadron analysis software

1. In your home directory, create an empty 'install' folder

2. Ensure both ROOTSYS and CLAS12ROOT are environment variables (see https://github.com/JeffersonLab/clas12root for details on installing CLAS12ROOT)

3. Run ./autogen.sh --prefix=/path/to/home/install

4. Run make

5. Run make install

6. Ensure both a lib/ and include/ directory has been created with dihadronana files stored inside
