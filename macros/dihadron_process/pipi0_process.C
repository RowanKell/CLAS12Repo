#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
//#include </work/clas12/users/gmat/install/include/dihadronana/SIDISKinematicsReco.h>
//R__LOAD_LIBRARY(libdihadronana.so)
#endif

//using namespace std;

int build_CLAS12SIDIS_A(
			const char * inputFile = "/work/clas12/users/gmat/dihadron/pipi0/data/raw/",
			const char * outputFile = "/work/clas12/users/gmat/dihadron/pipi0/data/raw/"
			)
{
  //---------------
  // Load libraries
  //---------------
  gSystem->Load("/work/clas12/users/gmat/install/lib/libdihadronana.so"); 
 
  //  SIDISKinematicsReco *ana = new SIDISKinematicsReco("");
 

  return 0;
}
