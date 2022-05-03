#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
R__LOAD_LIBRARY(libclas12ana.so)
#endif

using namespace std;

int pipi0_process(
			const char * inputFile = "/work/clas12/users/gmat/dihadron/pipi0/data/raw/",
			const char * outputFile = "/work/clas12/users/gmat/dihadron/pipi0/data/raw/"
			)
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("libclas12ana.so"); 
 
  //  SIDISKinematicsReco *ana = new SIDISKinematicsReco("");
 

  return 0;
}
