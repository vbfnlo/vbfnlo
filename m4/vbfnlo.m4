dnl most stuff in here is shamelessly adapted from Herwig++ or ThePEG

dnl some switches I think are useful

AC_DEFUN([VBFNLO_CHECK_KK],
[
dnl the check for the GSL library, which is needed for enable_kk=yes, is done later on in VBFNLO_CHECK_GSL
AC_MSG_CHECKING([whether KK should be included])
AC_ARG_ENABLE(kk,
        AC_HELP_STRING([--enable-kk],[--enable-kk to enable simulation of Kaluza-Klein resonances]),
        [],
        [enable_kk=no]
        )
AC_MSG_RESULT([$enable_kk])
AM_CONDITIONAL(WITH_KK,[test "x$enable_kk" = "xyes"])
])

AC_DEFUN([VBFNLO_CHECK_NLO],
[
AC_MSG_CHECKING([whether NLO should be enabled])
AC_ARG_ENABLE(NLO,
        AC_HELP_STRING([--enable-NLO],[--disable-NLO to disable next-to-leading order QCD]),
        [],
        [enable_NLO=yes]
        )
AC_MSG_RESULT([$enable_NLO])
AM_CONDITIONAL(WITH_NLO,[test "x$enable_NLO" = "xyes"])
])

AC_DEFUN([VBFNLO_CHECK_MPI],
[
AC_MSG_CHECKING([whether MPI should be enabled])
AC_ARG_ENABLE(MPI,
        AC_HELP_STRING([--enable-MPI],[--enable-MPI to use MPI parallelization]),
        [],
        [enable_MPI=no]
        )
AC_MSG_RESULT([$enable_MPI])
AM_CONDITIONAL(WITH_MPI,[test "x$enable_MPI" = "xyes"])
])


AC_DEFUN([VBFNLO_CHECK_QUAD],
[
AC_MSG_CHECKING([whether quad precision should be enabled])
AC_ARG_ENABLE(quad,
        AC_HELP_STRING([--enable-quad],[--enable-quad to enable quad precision for difficult phase space points]),
        [],
        [enable_quad=$supports_quad]
        )

dnl now let's check if quad precision works
  if (test "x$enable_quad" = "xyes"); then
dnl check whether we know how to compile in quad precision
    if (test "x$VBFFC" = "xgfortran"); then
      AM_FCFLAGS="$AM_FCFLAGS -DQREAL=REALPART -DQIMAG=IMAGPART -DQCMPLX=COMPLEX"
      QUADFCFLAGS="-fdefault-real-8 -DQUAD=1 -DKIND=2"
    elif (test "x$VBFFC" = "xifort"); then
      QUADFCFLAGS="-DDBLE=QEXT -DDIMAG=QIMAG -DDCONJG=QCONJG -DDCMPLX=QCMPLX -r16 -DQUAD=1 -DKIND=2"
    else 
      enable_quad=no
      AC_MSG_ERROR([unsupported quad precision with compiler ($VBFFC)])
    fi
    if (test "x$enable_quad" = "xyes"); then
      oldFCFLAGS="$FCFLAGS"
      FCFLAGS="$FCFLAGS $QUADFCFLAGS"
      AC_LANG_ASSERT(Fortran)
      AC_LINK_IFELSE(
        AC_LANG_SOURCE([[       double precision x
                                x=1e-20
                                x=1+x 
                                if (x.eq.1d0) then
                                  print*, "ERROR"
                                else
                                  print*, "OK"
                                endif
                                end
                       ]]),
        [if (./conftest$EXEEXT | grep -q "OK"); then
           AC_SUBST(QUADFCFLAGS)
           AC_MSG_RESULT([$enable_quad])
         else
           enable_quad=no
           AC_MSG_ERROR([quad precision not reached])
         fi ],
        [enable_quad=no
         AC_MSG_ERROR([compilation failed])]
      )
      FCFLAGS="$oldFCFLAGS"
    fi
  else
    AC_MSG_RESULT([$enable_quad])
  fi

AM_CONDITIONAL(WITH_QUAD,[test "x$enable_quad" = "xyes"])
])


dnl which processes to enable

AC_DEFUN([VBFNLO_ENABLE_PROCESSES],
[
AC_MSG_CHECKING([for processes to include])

AC_ARG_ENABLE(processes,
        AC_HELP_STRING([--enable-processes=LIST],[Comma-separated list of processes to enable. Options are: vbf, ggf, diboson, dibosonjet, triboson, tribosonjet, hjjj, qcdvvjj, qcdvjj, all_except_hexagons, all. Default is all_except_hexagons.]),
        [],
        [enable_processes=all_except_hexagons]
        )
AC_MSG_RESULT([$enable_processes])

oldIFS="$IFS"
IFS=","
for i in $enable_processes; do
    declare $i=yes
done
IFS="$oldIFS"

AM_CONDITIONAL(WITH_DIBOSON,[test "$diboson" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_GGF,[test "$ggf" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_QCDVV,[test "$qcdvvjj" -o "$all"])
AM_CONDITIONAL(WITH_DIBOSONJET,[test "$dibosonjet" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_TRIBOSON,[test "$triboson" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_QCDV,[test "$qcdvjj" -o "$qcdvvjj" -o "$all" ])
AM_CONDITIONAL(WITH_VBF,[test "$vbf" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_HJJJ,[test "$hjjj" -o "$all" -o "$all_except_hexagons"])
AM_CONDITIONAL(WITH_TRIBOSONJET,[test "$tribosonjet" -o "$all"])

dnl do we need HexLine?
AM_CONDITIONAL(WITH_HEXAGONS,[test "$tribosonjet" -o "$all" -o "$qcdvvjj"])
dnl do we need FermionLoops?
AM_CONDITIONAL(WITH_FERMIONLOOPS,[test "$diboson" -o "$dibosonjet" -o "$qcdvjj" -o "$qcdvvjj" -o "$all" -o "$all_except_hexagons"])
dnl do we need Pentagons with 4 partons; Hexagons are only needed for QCDVV
AM_CONDITIONAL(WITH_QCD_JJ_PLUSX,[test "$all" -o "$qcdvjj" -o "$qcdvvjj" ])
dnl check if quad precision is desirable
AM_COND_IF([WITH_FERMIONLOOPS],
  [supports_quad=yes],
  [
if [test "$all" -o "$qcdvjj" -o "$qcdvvjj" ]; then
  supports_quad=yes
else
  supports_quad=no
fi
  ])
])

dnl ##### ROOT #####
AC_DEFUN([VBFNLO_CHECK_ROOT],[

ROOTPATH=""
ROOTLIBPATH=""
ROOTLIBS=""
ROOTLDFLAGS=""
ROOTINCLUDE=""


AC_ARG_WITH(root,
        AC_HELP_STRING([--with-root=DIR],[location of ROOT installation]),
        [],
        [with_root=no])

if test "x$with_root" = "xno"; then
  AC_MSG_CHECKING([for ROOT])
  AC_MSG_RESULT([not required])
  HAS_ROOT="no"
else
  HAS_ROOT="yes"

  if test \( "x$with_root" = "xyes" \) -o \( "x$with_root" = "x" \); then
    with_root="$ROOTSYS"
  fi

  ROOTPATH="$with_root"

  if test -z "$ROOTPATH"; then
  	dnl assume root-config in system path
  	AC_PATH_PROG(ROOTCONFIG, root-config)
  else
  	AC_PATH_PROG(ROOTCONFIG, root-config, [""], ${ROOTPATH}/bin)
  fi

  AC_MSG_CHECKING([for ROOT])
  if test -z "${ROOTCONFIG}"; then
    if test -f $ROOTPATH/bin/root-config; then
        ROOTCONFIG=$ROOTPATH/bin/root-config
        AC_MSG_RESULT([$ROOTCONFIG])
    else
        AC_MSG_ERROR([--with-root was given, got $ROOTPATH but did not find $ROOTPATH/bin/root-config])
    fi
  else
    AC_MSG_RESULT([$ROOTCONFIG])
  fi

  AC_MSG_CHECKING([for ROOTLIBS])
  if test -z "$ROOTLIBS"; then
    ROOTLIBS=`$ROOTCONFIG --libs`
    ROOTLIBPATH=`$ROOTCONFIG --libdir`
    ROOTLDFLAGS="-L$ROOTLIBPATH -Wl,-rpath,$ROOTLIBPATH"
  fi
  AC_MSG_RESULT([$ROOTLIBS])

  AC_MSG_CHECKING([for ROOTINCLUDE])
  if test -z "$ROOTINCLUDE"; then
    ROOTINCLUDE="`$ROOTCONFIG --cflags`"
  fi
  AC_MSG_RESULT([$ROOTINCLUDE])

  oldLIBS="$LIBS"
  oldLDFLAGS="$LDFLAGS"
  oldCPPFLAGS="$CPPFLAGS"
  LIBS="$LIBS $ROOTLIBS"
  LDFLAGS="$LDFLAGS $ROOTLDFLAGS"
  CPPFLAGS="$CPPFLAGS $ROOTINCLUDE"

  AC_MSG_CHECKING([that ROOT works])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <TCanvas.h>]],
  	[[TCanvas c("c1", "");]])],[AC_MSG_RESULT([yes])],[AC_MSG_RESULT([no]) 
   	AC_MSG_ERROR([Use '--with-root=' to set the path to your ROOT installation.\
  If it doesn't work anyhow, you eventually have to set the ROOTSYS environment variable.])
  ])
  AC_LANG_POP([C++])

  LIBS="$oldLIBS"
  LDFLAGS="$oldLDFLAGS"
  CPPFLAGS="$oldCPPFLAGS"

  AC_SUBST(ROOTLIBS)
  AC_SUBST(ROOTLIBPATH)
  AC_SUBST(ROOTLDFLAGS)
  AC_SUBST(ROOTPATH)
  AC_SUBST(ROOTINCLUDE)

fi
AM_CONDITIONAL(WITH_ROOT, test x"$HAS_ROOT" = "xyes")
])

dnl ###### GSL ######
AC_DEFUN([VBFNLO_CHECK_GSL],
[
AC_MSG_CHECKING([for gsl location])
GSLINCLUDE=""
GSLLDFLAGS=""
GSLLIBS=""

AC_ARG_WITH(gsl,
        AC_HELP_STRING([--with-gsl=DIR],[location of gsl installation @<:@default=system libs@:>@]),
        [],
	  [with_gsl=system])

dnl check if GSL is needed
if test "x$enable_kk" = "xyes"; then

  if test "x$with_gsl" = "xyes"; then
    with_gsl="system"
  fi

  if test "x$with_gsl" = "xsystem"; then
      AC_LANG_PUSH([C])
	AC_MSG_RESULT([in system libraries])
	oldlibs="$LIBS"
	AC_SEARCH_LIBS(gsl_integration_qags,gsl,[],
			[
			AC_MSG_ERROR([Cannot find libgsl. Please install the GNU scientific library and header files into your system libraries or use --with-gsl=/path/to/GSL/ !])
			],
                  [-lgslcblas -lm]
		      )
	GSLLIBS="$LIBS -lgslcblas -lm "
	LIBS=$oldlibs
      AC_LANG_POP([C])
  elif test "x$with_gsl" = "xno"; then
	AC_MSG_RESULT([disabled])
	AC_MSG_ERROR([GSL libraries have been switched off, but they are needed for --enable_kk !])
  else
	if test "`uname -m`" = "x86_64" -a \( -e "$with_gsl/lib64/libgsl.a" -o -e "$with_gsl/lib64/libgsl.so" \) -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib64 -Wl,-rpath,$with_gsl/lib64 -lgslcblas -lgsl -lm"
            GSLLDFLAGS="-L$with_gsl/lib64 -Wl,-rpath,$with_gsl/lib64"
		GSLINCLUDE="-I$with_gsl/include"
	elif test \( -e "$with_gsl/lib/libgsl.a" -o -e "$with_gsl/lib/libgsl.so" \) -a -d "$with_gsl/include/gsl"; then
		AC_MSG_RESULT([found in $with_gsl])
		GSLLIBS="-L$with_gsl/lib -Wl,-rpath,$with_gsl/lib -lgslcblas -lgsl -lm"
            GSLLDFLAGS="-L$with_gsl/lib64 -Wl,-rpath,$with_gsl/lib64"
		GSLINCLUDE="-I$with_gsl/include"
	else
		AC_MSG_RESULT([not found])
		AC_MSG_ERROR([Can't find $with_gsl/lib/libgsl.a / $with_gsl/lib/libgsl.so or the headers in $with_gsl/include])
	fi
  fi
  GSLLIBS+=" -lstdc++"

else

  with_gsl="not required"
  AC_MSG_RESULT([not required])

fi

AM_CONDITIONAL([WITH_GSL], [test "x$with_gsl" != "xnot required"])


AC_SUBST(GSLINCLUDE)
AC_SUBST(GSLLDFLAGS)
AC_SUBST(GSLLIBS)
])

dnl #### LHAPDF ######
AC_DEFUN([VBFNLO_CHECK_LHAPDF],
[

AC_ARG_WITH(LHAPDF,
			AC_HELP_STRING([--with-LHAPDF=DIR],[location of LHAPDF installation]),
dnl			[  --with-LHAPDF=/path/to/LHAPDF/ or --without-LHAPDF to use internal PDF sets], 
			[],
			[with_LHAPDF=no])

if test "x$with_LHAPDF" = "xno"; then
	HAS_LHAPDF="no"
	LHAPDF_DIR=""
	LHAPDF_LDFLAGS=""
	LHAPDF_LIBS=""
	LHAPDF_PDFPATH=""
	AC_MSG_CHECKING([for LHAPDF])
	AC_MSG_RESULT([not required])
else
	HAS_LHAPDF="yes"
	if test "x$with_LHAPDF" = "xyes"; then
		LHAPDF_DIR=""
	else
		LHAPDF_DIR=$with_LHAPDF
	fi
  
	dnl search for the lhapdf-config script
	if test -z "$LHAPDF_DIR"; then
		dnl assume lhapdf-config in system path
		AC_PATH_PROG(lhaconfig, lhapdf-config)
	else
		AC_PATH_PROG(lhaconfig, lhapdf-config, [""], ${LHAPDF_DIR}/bin)
	fi

	if test -z "${lhaconfig}"; then
	    if test -z "$LHAPDF_DIR"; then
			AC_MSG_ERROR([--with-lhapdf was given, but can't find the "lhapdf-config" program in PATH]);
		else
			AC_MSG_ERROR([--with-lhapdf was given, but can't find the "lhapdf-config" program in ${LHAPDF_DIR}/bin]);
		fi
	else
	   dnl now see if LHAPDF is functional
       AC_MSG_CHECKING([for LHAPDF])
	   save_LDFLAGS="$LDFLAGS"
	   save_LIBS="$LIBS"

	   LDFLAGS="${LDFLAGS} -L`${lhaconfig} --libdir`"
	   LIBS="${LIBS} -lLHAPDF"
	   AC_LANG_PUSH(C++)
	   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern "C" { void initpdf_(int&); }]], 
						  [[int i = 1; initpdf_(i);]])], 
						  [lhaok='yes'], [lhaok='no'])
	   AC_LANG_POP()
	   LDFLAGS="$save_LDFLAGS"
	   LIBS="$save_LIBS"

	   if test "${lhaok}" = "yes"; then
		  LHAPDF_LDFLAGS="-L`${lhaconfig} --libdir`"
		  LHAPDF_LIBS="-lLHAPDF"
		  LHAPDF_PDFPATH="`${lhaconfig} --datadir`"
		  LHAPDF_VERSION="`${lhaconfig} --version | cut -d. -f1`"
		  LHAPDF_DIR="`${lhaconfig} --prefix`"
		  AC_MSG_RESULT([yes])
	   else
		  AC_MSG_RESULT([no])
		  AC_MSG_ERROR([--with-lhapdf was given, but can't find LHAPDF])
	   fi
	fi

	dnl output LHAPDF version
	AC_MSG_CHECKING([which LHAPDF version are we using])
	AC_MSG_RESULT([$LHAPDF_VERSION])
fi

AC_SUBST(LHAPDF_DIR)
AC_SUBST(LHAPDF_LIBS)
AC_SUBST(LHAPDF_LDFLAGS)
AC_SUBST(LHAPDF_PDFPATH)
AM_CONDITIONAL([WITH_LHAPDF], [test "x$HAS_LHAPDF" == "xyes"])
AM_CONDITIONAL([WITH_LHAPDF_6], [test "x$LHAPDF_VERSION" == "x6"])
AC_DEFINE_UNQUOTED([LHAPDFPATH],["$LHAPDF_PDFPATH"],[PDF directory])
])

dnl #### HEPMC ######
AC_DEFUN([VBFNLO_CHECK_HEPMC],
[
AC_MSG_CHECKING([for HepMC])
AC_ARG_WITH(hepmc,
        AC_HELP_STRING([--with-hepmc=DIR],[location of HepMC installation]),
        [],
        [with_hepmc=no])

if test "x$with_hepmc" = "xno"; then
  AC_MSG_RESULT([not required])
else

HAS_HEPMC="yes"

if test "x$with_hepmc" != "xyes"; then
	if test "`uname -m`" = "x86_64" -a -d "$with_hepmc/lib64" ; then
		HEPMCLIBS="-L$with_hepmc/lib64 -Wl,-rpath,$with_hepmc/lib64 -lHepMC"
		HEPMCLDFLAGS="-L$with_hepmc/lib64 -Wl,-rpath,$with_hepmc/lib64"
	else
		HEPMCLIBS="-L$with_hepmc/lib -Wl,-rpath,$with_hepmc/lib -lHepMC"
		HEPMCLDFLAGS="-L$with_hepmc/lib -Wl,-rpath,$with_hepmc/lib"
	fi
        HEPMC_DIR="$with_hepmc"
	HEPMCINCLUDE="-I$with_hepmc/include"
else
	HEPMCLIBS="-lHepMC"
	HEPMCLDFLAGS=""
	HEPMC_DIR=""
	HEPMCINCLUDE=""
fi
AC_LANG_PUSH([C++])
oldLIB="$LIBS"
oldLDFLAGS="$LDFLAGS"
oldCPPFLAGS="$CPPFLAGS"
LIBS="$LIBS $HEPMCLIBS"
LDFLAGS="$LDFLAGS $HEPMCLDFLAGS"
CPPFLAGS="$CPPFLAGS $HEPMCINCLUDE"
AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([[#include <HepMC/GenEvent.h>]],
    	            [[HepMC::write_HepMC_IO_block_begin(std::cout);]]
                   )],
   [AC_MSG_RESULT([yes])],
   [AC_MSG_RESULT([no]) 
    AC_MSG_ERROR([Cannot find libHepMC. Please install the HepMC library and header files])]
)
LIBS=$oldLIB
LDFLAGS=$oldLDFLAGS
CPPFLAGS=$oldCPPFLAGS
AC_LANG_POP([C++])
fi
AC_SUBST(HEPMC_DIR)
AC_SUBST(HEPMCLIBS)
AC_SUBST(HEPMCLDFLAGS)
AC_SUBST(HEPMCINCLUDE)
AM_CONDITIONAL([WITH_HEPMC], [test "x$HAS_HEPMC" == "xyes"])
])


dnl #### LOOPTOOLS ######

AC_DEFUN([VBFNLO_CHECK_LOOPTOOLS],
[
AC_MSG_CHECKING([for LOOPTOOLS])

HAS_LOOPTOOLS="yes"
LOOPTOOLS_DIR=""
LOOPTOOLS_LIBDIR=""
AC_ARG_WITH(LOOPTOOLS,
            [  --with-LOOPTOOLS=/path/to/LOOPTOOLS/ to enable the calculation of electroweak corrections], 
            [if test -n "$with_LOOPTOOLS" -a "x$with_LOOPTOOLS" != "xyes" -a "x$with_LOOPTOOLS" != "xno"; then 
               LOOPTOOLS_DIR="$with_LOOPTOOLS"; 
             elif test "x$with_LOOPTOOLS" = "xno"; then 
               HAS_LOOPTOOLS="no"; 
             fi])


if test -n "$LOOPTOOLS_DIR"; then
    HAS_LOOPTOOLS="yes"

  if test -x "$LOOPTOOLS_DIR/lib"; then
    LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib
  elif test -x "$LOOPTOOLS_DIR/lib64"; then
    LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib64
  elif test -x "$LOOPTOOLS_DIR/$(tcsh -cf 'echo $HOSTTYPE')"; then
    LOOPTOOLS_DIR=$LOOPTOOLS_DIR/$(tcsh -cf 'echo $HOSTTYPE')
    if test -x "$LOOPTOOLS_DIR/lib"; then
      LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib
    elif test -x "$LOOPTOOLS_DIR/lib64"; then
      LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib64
    fi
  elif test -x "$LOOPTOOLS_DIR/$(tcsh -cf 'uname -m')-$(tcsh -cf 'uname -s')"; then
    LOOPTOOLS_DIR=$LOOPTOOLS_DIR/$(tcsh -cf 'uname -m')-$(tcsh -cf 'uname -s')
    if test -x "$LOOPTOOLS_DIR/lib"; then
      LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib
    elif test -x "$LOOPTOOLS_DIR/lib64"; then
      LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR/lib64
    fi
  else
    LOOPTOOLS_LIBDIR=$LOOPTOOLS_DIR
    LOOPTOOLS_DIR=$(echo "$LOOPTOOLS_LIBDIR" | sed 's/lib64//g')
    LOOPTOOLS_DIR=$(echo "$LOOPTOOLS_DIR" | sed 's/lib//g')
  fi

  LOOPTOOLS_LDFLAGS="-L$LOOPTOOLS_LIBDIR"
  LOOPTOOLS_LIBS="-looptools"

 AC_MSG_CHECKING([link to LoopTools])
 oldlibs="$LIBS"
 LIBS="$LOOPTOOLS_LDFLAGS $LOOPTOOLS_LIBS"
 AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([],
                    [[        call ffini()
                                test = A0(100D0)]])], 
   [LTVERSION="OLD"],
 [
   AC_LINK_IFELSE(
     [AC_LANG_PROGRAM([],
                      [[        call ltini()
                                test = A0(100D0)]])], 
     [LTVERSION="NEW"],
  AC_MSG_ERROR(["Could not link with LoopTools. Are you sure the path given is correct and that the same compiler (gfortran/ifort) is used for VBFNLO and LoopTools?"]))
  ]
 )

 LIBS=$oldlibs

 oldLTDIR=$LOOPTOOLS_DIR
 LOOPTOOLS_DIR=$(echo "$LOOPTOOLS_DIR" | sed "s/$(tcsh -cf 'uname -m')-$(tcsh -cf 'uname -s')//g")
 LOOPTOOLS_DIR=$(echo "$LOOPTOOLS_DIR" | sed "s/$(tcsh -cf 'echo $HOSTTYPE')//g")
 if test -n "$(fgrep 'ff2c' $LOOPTOOLS_DIR/makefile)"; then
    AM_FCFLAGS="$AM_FCFLAGS -ff2c"
 fi
 LOOPTOOLS_DIR=$oldLTDIR


  { $as_echo "$as_me:${as_lineno-$LINENO}: result: yes" >&5
    $as_echo "okay" >&6; }

 if test "${enable_shared}" != "no"; then
  if objdump --reloc $LOOPTOOLS_LIBDIR/libooptools.a | grep -q R_X86_64_32S; then
   AC_MSG_ERROR(["The version of LoopTools you have linked to is statically linked.  Please re-run configure with the option '--enable-shared=no'."])
  fi
 fi

else
    HAS_LOOPTOOLS="no"
  { $as_echo "$as_me:${as_lineno-$LINENO}: result: not required" >&5
    $as_echo "not required" >&6; }
fi


AM_CONDITIONAL([NEW_LT], [test "$LTVERSION" == "NEW"])
AM_CONDITIONAL([OLD_LT], [test "$LTVERSION" == "OLD"])
AC_SUBST(LOOPTOOLS_DIR)
AC_SUBST(LOOPTOOLS_LIBS)
AC_SUBST(LOOPTOOLS_LDFLAGS)
AM_CONDITIONAL([WITH_LOOPTOOLS], [test "$HAS_LOOPTOOLS" == "yes"])
])


dnl #### FEYNHIGGS ######

AC_DEFUN([VBFNLO_CHECK_FEYNHIGGS],
[
AC_MSG_CHECKING([for FEYNHIGGS])

HAS_FEYNHIGGS="yes"
FEYNHIGGS_DIR=""
FEYNHIGGS_LIBDIR=""
AC_ARG_WITH(FEYNHIGGS,[  --with-FEYNHIGGS=/path/to/FEYNHIGGS/ to enable calculation of MSSM Higgs sector by FeynHiggs], [if test -n "$with_FEYNHIGGS" -a "x$with_FEYNHIGGS" != "xyes" -a "x$with_FEYNHIGGS" != "xno"; then FEYNHIGGS_DIR="$with_FEYNHIGGS"; elif test "x$with_FEYNHIGGS" = "xno"; then HAS_FEYNHIGGS="no"; fi])


if test -n "$FEYNHIGGS_DIR"; then

  if test -x "$FEYNHIGGS_DIR/lib"; then
    FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib
  elif test -x "$FEYNHIGGS_DIR/lib64"; then
    FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib64
  elif test -x "$FEYNHIGGS_DIR/$(tcsh -cf 'echo $HOSTTYPE')"; then
    FEYNHIGGS_DIR=$FEYNHIGGS_DIR/$(tcsh -cf 'echo $HOSTTYPE')
    if test -x "$FEYNHIGGS_DIR/lib"; then
      FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib
    elif test -x "$FEYNHIGGS_DIR/lib64"; then
      FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib64
    fi
  elif test -x "$FEYNHIGGS_DIR/$(tcsh -cf 'uname -m')-$(tcsh -cf 'uname -s')"; then
    FEYNHIGGS_DIR=$FEYNHIGGS_DIR/$(tcsh -cf 'uname -m')-$(tcsh -cf 'uname -s')
    if test -x "$FEYNHIGGS_DIR/lib"; then
      FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib
    elif test -x "$FEYNHIGGS_DIR/lib64"; then
      FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR/lib64
    fi
  else
    FEYNHIGGS_LIBDIR=$FEYNHIGGS_DIR 
    FEYNHIGGS_DIR=$(echo "$FEYNHIGGS_LIBDIR" | sed 's/lib64//g')
    FEYNHIGGS_DIR=$(echo "$FEYNHIGGS_DIR" | sed 's/lib//g')
  fi

 FEYNHIGGS_LDFLAGS="-L$FEYNHIGGS_LIBDIR"
 FEYNHIGGS_LIBS="-lFH"

 AC_MSG_CHECKING([link to FeynHiggs])
 oldlibs="$LIBS"
 LIBS="$FEYNHIGGS_LDFLAGS $FEYNHIGGS_LIBS"
AC_LINK_IFELSE(
   [AC_LANG_PROGRAM([],
                    [[           call FHSetFlags(0,1,1,1,1,1,1,1,1,1)]])], 
   [ HAS_FEYNHIGGS="yes"],
 [ AC_MSG_ERROR(["Could not link with FeynHiggs.  Are you sure the path given is correct and that the same compiler (gfortran/ifort) is used for VBFNLO and FeynHiggs?"])
  ]
 )
 LIBS=$oldlibs


  { $as_echo "$as_me:${as_lineno-$LINENO}: result: yes" >&5
    $as_echo "okay" >&6; }

 FHSTRING=$( LC_ALL=C grep --binary-files=text -oP 'FeynHiggs \d+.\d+.\d+' $FEYNHIGGS_LIBDIR/libFH.a  | sed -e 's/FeynHiggs //' )
 FHVERSION=$( echo "$FHSTRING" | $AWK -F'.' '{print $[1]*10000+$[2]*100+$[3]}' )
 dnl reconstruct string for verification
 FHSTRING="$((FHVERSION/10000)).$((FHVERSION/100-FHVERSION/10000*100)).$((FHVERSION-FHVERSION/100*100))"

 if test "x$FHSTRING" = "x"; then
  AC_MSG_ERROR(["Sorry, we don't recognise the version of FeynHiggs you are using."])
 fi
 if test $FHVERSION -lt 20600; then
  AC_MSG_ERROR(["Sorry, your version of FeynHiggs is too old ($FHSTRING, at least 2.6.0 required)."])
 fi

 if test "${enable_shared}" != "no"; then
  if objdump --reloc $FEYNHIGGS_LIBDIR/libFH.a | grep -q R_X86_64_32S; then
   AC_MSG_ERROR(["The version of FeynHiggs you have linked to is statically linked.  Please re-run configure with the option '--enable-shared=no'."])
  fi
 fi

else
    HAS_FEYNHIGGS="no"
  { $as_echo "$as_me:${as_lineno-$LINENO}: result: not required" >&5
    $as_echo "not required" >&6; }
fi


AC_SUBST(FHVERSION)
AC_SUBST(FEYNHIGGS_DIR)
AC_SUBST(FEYNHIGGS_LIBS)
AC_SUBST(FEYNHIGGS_LDFLAGS)
AM_CONDITIONAL([WITH_FEYNHIGGS], [test "x$HAS_FEYNHIGGS" == "xyes"])
])


dnl #### check, which fortran compiler we are using #####
AC_DEFUN([VBFNLO_CHECK_FC],
[
  AC_MSG_CHECKING([which Fortran compiler we are using])

  dnl TODO: rewrite with the right autotools functions
  check_g77=`echo $FC | grep g77`
  check_gfortran=`echo $FC | grep gfortran || echo $FC | grep mpifort`
  check_ifort=`echo $FC | grep ifort`

  if test -n "$check_g77" ; then
     AC_MSG_ERROR([g77 is no longer supported by VBFNLO. Please use gfortran or ifort!])
  elif test -n "$check_gfortran" ; then
	 GFORTRAN_VERSION="`$FC -dumpversion`"
	 AS_VERSION_COMPARE([${GFORTRAN_VERSION}], [4.2], 
		 AC_MSG_ERROR("You need gfortran at version at least 4.2 to compile VBFNLO. Please update gfortran or use ifort.")
	 , [], [])

     VBFFC="gfortran"
     AM_FCFLAGS="$AM_FCFLAGS -fno-automatic -ffixed-line-length-none -ffree-line-length-none -fimplicit-none -fbacktrace"
     GFORTRAN_48_FIX_LOOPS="-fautomatic"
     GFORTRAN_48_FIX_GGF="-O1"
  elif test -n "$check_ifort" ; then
     VBFFC="ifort"
     dnl for ifort the line length for free-format source files (.F90/.f90) is per default unlimited.
     dnl -extend_source allows line lengths in fixed-form files (.F/.f) of 132 chars.
     dnl in ifort: max. 511 continuation lines for fixed-form, 255 continuation lines for free-form
     AM_FCFLAGS="$AM_FCFLAGS -save -extend_source -traceback"
     GFORTRAN_48_FIX_LOOPS="-auto"
     GFORTRAN_48_FIX_GGF="-O1"
  else
     VBFFC="unrecognized, assuming gfortran behaviour"
     AM_FCFLAGS="$AM_FCFLAGS -fno-automatic -ffixed-line-length-none -ffree-line-length-none -fimplicit-none -fbacktrace"

     GFORTRAN_48_FIX_LOOPS="-fautomatic"
     GFORTRAN_48_FIX_GGF="-O1"
  fi

dnl gfortran 4.8 has a memory and compile time problem with vbfnlo in the loop routines.
dnl In order to circumvent this some loop libraries are compiled without an implicit "save" on variables, 
dnl which is fine in these routines.
dnl Additionally, the compile time for some GGF routines increase roughly by a factor of 60 with -O2,
dnl so they are compiled with -O1 (which at least improves the time by a factor of 10).
  AC_SUBST(GFORTRAN_48_FIX_LOOPS)
  AC_SUBST(GFORTRAN_48_FIX_GGF)

  AC_MSG_RESULT([$VBFFC])

  AC_FC_LIBRARY_LDFLAGS

])

dnl check for MadGraph comparison mode and set flags
AC_DEFUN([VBFNLO_CHECK_MADGRAPH],
[
AC_MSG_CHECKING([whether MadGraph comparisons should be included])
AC_ARG_ENABLE(madgraph,
        AC_HELP_STRING([--enable-madgraph],[--enable-madgraph to include code for MadGraph comparisons]),
        [],
        [enable_madgraph=no]
        )
AC_MSG_RESULT([$enable_madgraph])
AM_CONDITIONAL(WITH_MADGRAPH,[test "x$enable_madgraph" = "xyes"])
])

dnl check for git checkout

AC_DEFUN([VBFNLO_VERSIONSTRING],
[
AC_CHECK_PROG(have_git,[git],[yes],[no])
AC_MSG_CHECKING([whether a git version string is available])
if test -e $srcdir/.git; then
    if test x"$have_git" != x"yes" ; then
        AC_MSG_ERROR([Need git installed if working with a git checkout (.git is present).])
    fi
else
    have_git="no"
fi
AC_MSG_RESULT([$have_git])
])

dnl set naming conventions for Fortran-C linking

AC_DEFUN([BEFORE_CFORTRANLINK],
[
  oldfcflags="$FCFLAGS"
  FCFLAGS="$AM_FCFLAGS"
])
AC_DEFUN([VBFNLO_CFORTRANLINK],
[
  AC_REQUIRE([BEFORE_CFORTRANLINK])
  AC_FC_WRAPPERS
  FCFLAGS=$oldfcflags
])

dnl check for debugging mode and set compiler flags

AC_DEFUN([VBFNLO_SET_FLAGS],
[

  if (test "x$VBFFC" = "xgfortran"); then
    AM_FCFLAGS="$AM_FCFLAGS -J\$(top_builddir)/include"
  elif (test "x$VBFFC" = "xifort"); then
    AM_FCFLAGS="$AM_FCFLAGS -Wc,-module,\$(top_builddir)/include"
  fi

AC_MSG_CHECKING([for debugging mode])
AC_ARG_ENABLE(debug,
        AC_HELP_STRING([--enable-debug],[--enable-debug to enable more compile warnings]),
        [],
        [enable_debug=no]
        )
AC_MSG_RESULT([$enable_debug])
AM_CONDITIONAL(VBFNLO_DEBUG,[test "x$enable_debug" = "xyes"])


if test "x$enable_debug" == "xyes" ; then
  if (test "x$VBFFC" = "xgfortran"); then
    AM_FCFLAGS="$AM_FCFLAGS -Wall -Wextra"
  elif (test "x$VBFFC" = "xifort"); then
    AM_FCFLAGS="$AM_FCFLAGS -warn"
  fi
  AM_CFLAGS="$AM_CFLAGS -Wall -Wextra"
  AM_CXXFLAGS="$AM_CXXFLAGS -Wall -Wextra"
fi
  AM_FCFLAGS="$AM_FCFLAGS -g -I\$(top_builddir)/include -O2"
  AM_CPPFLAGS="$AM_CPPFLAGS -g -I\$(top_builddir)/include"
  AM_CXXFLAGS="$AM_CXXFLAGS -g"
  AM_LDFLAGS="$AM_LDFLAGS $LDFLAGS $ROOTLDFLAGS $GSLLDFLAGS $LHAPDF_LDFLAGS $HEPMCLDFLAGS $LOOPTOOLS_LDFLAGS $FEYNHIGGS_LDFLAGS" 

  AC_SUBST(AM_CPPFLAGS)
  AC_SUBST(AM_CFLAGS)
  AC_SUBST(AM_CXXFLAGS)
  AC_SUBST(AM_FCFLAGS)
  AC_SUBST(AM_LDFLAGS)

  AC_SUBST([F77],[$FC])
  AC_SUBST([FFLAGS],[$FCFLAGS])
  AC_SUBST([AM_FFLAGS],[$AM_FCFLAGS])

  AM_CONDITIONAL(WITH_LOOPS,[test "x$enable_NLO" == "xyes" -o "$ggf" -o "$diboson" -o "$all" -o "$all_except_hexagons"])
])


dnl check for presence of pdflatex for Manual

AC_DEFUN([VBFNLO_CHECK_LATEX],
[
  AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
  if test -z "$PDFLATEX"; then
    AC_MSG_WARN([Unable to create PDF version of the user manual.])
  fi
  AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")
])


AC_DEFUN([VBFNLO_SUMMARY],
[
if test "x$with_LHAPDF" == "x" ; then
  with_LHAPDF="no"
fi
if test "x$with_LOOPTOOLS" == "x" ; then
  with_LOOPTOOLS="no"
fi
if test "x$with_FEYNHIGGS" == "x" ; then
  with_FEYNHIGGSversion="no"
else
  with_FEYNHIGGSversion="$with_FEYNHIGGS (version $FHSTRING)"
fi
if test "x$with_hepmc" == "x" ; then
  with_hepmc="no"
fi

echo "=================================================="
echo ""
echo "  $PACKAGE_STRING configuration"
echo ""
echo "  Processes included:"
echo "  $enable_processes"
echo ""
echo "  NLO QCD           $enable_NLO"
echo "  KK resonances     $enable_kk"
echo "  Compile warnings  $enable_debug"
echo "  Quad precision    $enable_quad"
if test "x$enable_madgraph" = "xyes" ; then
echo "  MadGraph code     included"
fi
echo "  MPI               $enable_MPI"
echo ""
echo "  LHAPDF            $with_LHAPDF"
echo "  LOOPTOOLS         $with_LOOPTOOLS"
echo "  FEYNHIGGS         $with_FEYNHIGGSversion"
echo "  ROOT              $with_root"
echo "  GSL               $with_gsl"
echo "  HEPMC             $with_hepmc"

echo ""
echo "  FORTRAN compiler  $VBFFC"

echo ""
echo "=================================================="
])

