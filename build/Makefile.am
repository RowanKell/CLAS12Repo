AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(CLAS12ROOT)/hipo4 \
  -I$(CLAS12ROOT)/Clas12Banks \
  -I$(CLAS12ROOT)/Clas12Root \
  -I$(ROOTSYS)/include \
  -I$(includedir) 

lib_LTLIBRARIES = \
   libclas12ana.la

AM_CXXFLAGS = -Wall -Werror -Wno-int-in-bool-context -msse2 -std=c++14

AM_LDFLAGS = \
  -L$(CLAS12ROOT)/lib \
  -L$(HIPO)/lib \
  -L$(libdir)

libclas12ana_la_LIBADD = \
  -lClas12Banks \
  -lClas12Root \
  -lHipo4 \
  -llz4

pkginclude_HEADERS = \
  ../src/Kinematics.h \
  ../src/Settings.h \
  ../src/SIDISKinematicsReco.h \
  ../src/SIDISParticle.h \
  ../src/Constants.h \
  ../src/PostProcess.h \
  ../src/SIDISParticlev1.h \
  ../src/FiducialCuts.h	

libclas12ana_la_SOURCES = \
  ../src/Kinematics.C \
  ../src/Settings.C \
  ../src/SIDISKinematicsReco.C \
  ../src/SIDISParticle.C \
  ../src/PostProcess.C \
  ../src/SIDISParticlev1.C \
  ../src/FiducialCuts.C

%_Dict.C: ../src/%.h ../src/%LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(AM_CPPFLAGS) $^
#	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^	

%_Dict_rdict.pcm: ../src/%_Dict.C ;

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.c

testexternals_LDADD = \
  libclas12ana.la

testexternals.c:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
