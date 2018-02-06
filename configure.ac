dnl name / version / bugreports / tarball name
AC_INIT([VBFNLO],[3.0.0beta5],[vbfnlo@itp.kit.edu],[VBFNLO])

dnl sanity check
AC_CONFIG_SRCDIR([src/vbfnlo_main.F90])

AC_CONFIG_MACRO_DIR([m4])

dnl place build and install tools in here
AC_CONFIG_AUX_DIR([config])

dnl fix warning of automake 1.12
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

dnl check for external compiler flags
if test "x$FCFLAGS" = "x" ; then
  dnl here one could add compiler flags which are appended after AM_FCFLAGS or library-specific FFLAGS
  dnl => overrides them, as external FC=xxx would do
  FCFLAGS=""
fi
FFLAGS=$FCFLAGS
if test "x$CXXFLAGS" = "x" ; then
    CXXFLAGS="-O2"
fi
if test "x$CCFLAGS" = "x" ; then
    CCFLAGS="-O2"
fi

AC_CHECK_TOOL(SWVERS, sw_vers)
if test x$SWVERS != x; then
  PROD_NAME=$($SWVERS -productName | cut -f 2 -d:)
fi
AM_CONDITIONAL(WITH_OSX, [test "$PROD_NAME" == "Mac OS X"])
if test "$PROD_NAME" == "Mac OS X"; then
  MACOSX_DEPLOYMENT_TARGET=$($SWVERS -productVersion | cut -f 1,2 -d.)
  AC_MSG_NOTICE([MACOSX_DEPLOYMENT_TARGET = $MACOSX_DEPLOYMENT_TARGET])
  AM_LDFLAGS="$AM_LDFLAGS -Wl,-undefined,dynamic_lookup"
  AM_CXXFLAGS="$AM_CXXFLAGS -Dunix"
  dnl fix tons of warnings for current gfortran
  if test "$MACOSX_DEPLOYMENT_TARGET" == "10.5" ; then
    AM_FCFLAGS="$AM_FCFLAGS -gdwarf-2"     
  fi
fi

dnl initialization
AC_LANG([Fortran])
AM_INIT_AUTOMAKE([1.9 gnu -Wall -Wno-portability subdir-objects foreign])
m4_ifdef([AM_EXTRA_RECURSIVE_TARGETS],
	 [AM_EXTRA_RECURSIVE_TARGETS([dependencylist])
          VBFNLO_DEPENDENCY_TRACKING()
	 ]
	)
AC_PROG_LIBTOOL([disable-static])
AC_PROG_F77()
AC_PROG_FC([], 95)
if test "x$FC" == "x" ; then
  AC_MSG_ERROR([No Fortran compiler found.])
fi
AC_PROG_INSTALL

dnl for root interfaces and kk utilities
AC_PROG_CXX

dnl to setup include dir
AC_PROG_LN_S

dnl check which fortran compiler is used and set compiler specific flags
VBFNLO_CHECK_FC


dnl check which processes to compile

VBFNLO_ENABLE_PROCESSES


dnl check for libraries

VBFNLO_CHECK_LHAPDF

VBFNLO_CHECK_LOOPTOOLS

VBFNLO_CHECK_FEYNHIGGS

VBFNLO_CHECK_ROOT

VBFNLO_CHECK_HEPMC


dnl check for requested VBFNLO features

VBFNLO_CHECK_KK

VBFNLO_CHECK_NLO

VBFNLO_CHECK_MPI

VBFNLO_CHECK_MADGRAPH


dnl check if quad precision works 
dnl (only tested if explicitly requested or if hexagon processes are included, see VBFNLO_ENABLE_PROCESSES)

VBFNLO_CHECK_QUAD


dnl the check for GSL depends on VBFNLO_CHECK_KK
VBFNLO_CHECK_GSL


dnl check if this is a developer version (in this case display GIT revision)
VBFNLO_VERSIONSTRING

dnl set naming conventions for Fortran-C linking
VBFNLO_CFORTRANLINK

dnl set compiler flags, depending on requested options
VBFNLO_SET_FLAGS

dnl check for presence of pdflatex for Manual
VBFNLO_CHECK_LATEX


dnl Makefiles to be generated
AC_CONFIG_FILES([include/Makefile
	         helas/Makefile
		 utilities/Makefile
		 PDFsets/Makefile
		 loops/Makefile
		 loops/Utilities/Makefile
		 loops/ggf_loops/Makefile
		 loops/Diagrams/Boxline/Makefile
		 loops/Diagrams/FermionLoops/Makefile		 
		 loops/Diagrams/HexLine/Makefile
		 loops/Diagrams/Penline/Makefile
		 loops/Diagrams/PenBox/Makefile
		 loops/Diagrams/HexBox/aa/Dire/Makefile
		 loops/Diagrams/HexBox/aa/Cross/Makefile
		 loops/Diagrams/HexBox/aq/Dire/Makefile
		 loops/Diagrams/HexBox/aq/Cross/Makefile
		 loops/Integrals/Makefile
		 loops/Integrals/FF/Makefile
		 loops/TenRed/B/Makefile	
		 loops/TenRed/C/Makefile
		 loops/TenRed/D/Makefile
		 loops/TenRed/E/Makefile
		 loops/TenRed/F/Makefile
		 phasespace/Makefile
		 amplitudes/Makefile
		 amplitudes/diboson/Makefile
		 amplitudes/triboson/Makefile
		 amplitudes/Spin2/Makefile
		 amplitudes/ggf/Makefile
		 amplitudes/qcdwz/Makefile
		 amplitudes/vvj/Makefile	
		 amplitudes/vvjj/Makefile
		 amplitudes/hjjj/Makefile	
         regress/Makefile
		 lib/Makefile
		 src/Makefile
                 doc/Makefile 
		 Makefile])

AC_CONFIG_HEADERS([utilities/VBFNLOConfig.h])

AC_OUTPUT

VBFNLO_SUMMARY
