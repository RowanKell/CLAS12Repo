# CLAS12Analysis
---
Analysis tool for reading hipo4 files with clas12root
# Setup Tutorial
---
1. Ensure both **ROOTSYS** and **CLAS12ROOT** are environment variables (see https://github.com/JeffersonLab/clas12root for details on installing CLAS12ROOT). For ROOT, we need 6.24.06, which should automatically install with 'module load clas12/pro'

2. Load in the CLAS12Analysis git
        git clone https://www.github.com/Gregtom3/CLAS12Analysis


3. Enter into the build directory

        cd build
        

4. Run 

        ./autogen.sh --prefix=/path/to/CLAS12Analysis/


5. Run `make` within the build directory

6. Run `make install` within the build directory

7. Leave the build directory and check to see that a lib and include directory have been created.

8. In your home directory, edit the .cshrc, add the following line

        set LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/CLAS12Analysis/lib
        
9. Reopen your terminal or `source ~/.cshrc` to have the LD_LIBRARY_PATH appended to


10. Keep in mind that, as of the current git build, I have not figured out how to `#include "<INSERT-CLAS12ANA-HEADER.h"` in programs such as **macros/dihadron_process/pipi0_process.C**. You will need to edit the path to the corresponding header files yourself. This also goes for slurmJobs in **slurm/**

# Structure
---
- **build/**
    - Directory containing necessary files for building/making the clas12ana.so library
- **data/**
    - Directory for storing ".root" files after processing/post-processing hipo4 files. The ".root" files stored here typically contain TTrees of particle kinematics, event variables, and other postprocessing features (such as diphoton mass for pi0 decays)
- **include/**
    - Directory created after running 'make' in "build". See **Setup** for instructions
- **lib/**
    - Same as above
- **macros/**
    - Directory which hosts two types of analysis scripts. All of these scripts are intended for the user to edit and run themselves. There is no other directory (other than **src/**) which the user must edit files within to conduct their desired analysis.
        - *Processing Scripts*: These C++ macros are run via [clas12root](https://github.com/JeffersonLab/clas12root) (ex: `clas12root pipi0_process.root`). The macro creates a `SIDISKinematicsReco` object from the compiled clas12ana.so library, and serves as the backbone of the analysis. A `Settings` object is additionally created, allowing the user to set cuts for the analysis using built-in functions. This object is also used to specify which banks the user wants to analyze (Monte Carlo and/or Reco). Lastly, after all events within the hipo4 file had been analyzed, the user can set a post processing method. This will trigger the `SIDISKinematicsReco` object to reload the reconstructed event TTree to perform additional analysis, storing the results in *tree_postprocess*. Some scripts are listed with the tag `postprocess_only`. These scripts will quickly load up a pre-existing reconstructed event TTree to run a newly defined post-processing method. Post-processing methods can be found in *src/PostProcess.C*
        - *Analysis Scripts*: These C++ macros are (typically) *.ipynb* files. They read in the TTree's (*tree_MC* , *tree_reco*, *tree_postprocess* ) created by the aforementioned processing scripts. 
- **slurm/**
    - Directory containing the shell script `./makeJobs.sh`. This shell script is edited by the user to run parallel Processing Scripts in **macros/**. The script is told which directory contains the .hipo files of interest, then loops over each file in the list, creating a slurm job for each. In the end `./makeJobs.sh` creates a new `./runJobs.sh` file, which, when executed, sends the slurm jobs to the computing cluster
- **src/**
    - Directory containing several .C and .h files for conducting the analysis. Workflow follows `SIDISKinematicsReco.C`, which loops over events in the `process_events` method. Any user edits to this file will not take into effect until the library is recompiled with make
- **tutorials/**
    - Directory containing some basic tutorial scripts for processing hipo files using the clas12ana.so library.
---
# Tutorials

### 1. Processing a small hipo file


To get started using the library, `cd tutorials/` and open up the *dihadron_tutorial.C* file with a text editor. Here, you can specify the hipo4 file desired for the analysis (its already setup for you at this stage), the cuts you want to apply, and the postprocessing method you want to complete. The file has plenty of comments to get you familiarized with how the cuts are set. Additional methods do exist in *src/Settings.h*, so feel free to check those out and play around with them.

To run the tutorial script, type...

        clas12root dihadron_tutorial.C
        
While running, the program uses functions native to `clas12root` to read hipo banks event-by-event. The program identifies the scattered electron to skip events which do not satisfy $Q^{2}$, $y$, and $W$ cuts. Additional particle specific cuts are made as well. So long as an event passes the reconstructed particle/event cuts, analysis from the monte carlo banks are also saved.

Using root, we can then open the `output.root` file to check out the resulting TTrees. I'd recommend printing each of the trees to see what exactly was stored. More information on how `tree_reco` and `tree_MC` store their data can be found by digging through **src/SIDISKinematicsReco.C**.

        root output.root
        tree_MC->Print()
        tree_reco->Print()
        tree_postprocess->Print()
        
If you have graphics enabled, you can plot the diphoton mass stored in `tree_postprocess`. The parameter (flag) is defined in  

        tree_postprocess->Draw("Mgg","flag==1")

---
Author: Gregory Matousek

Contact: gregory.matousek@duke.edu
