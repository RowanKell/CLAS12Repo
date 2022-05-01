AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(ROOTSYS)/include \
  -I$(CLAS12ROOT)/hipo4 \
  -I$(CLAS12ROOT)/Clas12Banks \
  -I$(CLAS12ROOT)/Clas12Root

lib_LTLIBRARIES = \
   libdihadronana.la

AM_CXXFLAGS = -Wall -Werror -msse2

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(CLAS12ROOT)/lib \
  -L$(HIPO)/lib

libdihadronana_la_LIBADD = \
  -lClas12Banks \
  -lClas12Root \
  -lHipo4 \
  -llz4

pkginclude_HEADERS = \
  Kinematics.h \
  Settings.h \
  SIDISKinematicsReco.h

libdihadronana_la_SOURCES = \
  Kinematics.C \
  Settings.C \
  SIDISKinematicsReco.C

%_Dict.C: %.h %LinkDef.h
	rootcint -f -std=c++11 $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.C ;

################################################
# linking tests

noinst_PROGRAMS = testexternals

BUILT_SOURCES = \
  testexternals.c

testexternals_LDADD = \
  libdihadronana.la

testexternals.c:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
