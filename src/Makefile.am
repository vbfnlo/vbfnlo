ADDSOURCES =
ADDFFLAGS = 
ADDLDADD =

if WITH_LOOPS
ADDSOURCES += FFini.F90
if WITH_QUAD
ADDFFLAGS += -DWITH_QUAD
endif
endif

if WITH_NLO
ADDFFLAGS += -DWITH_NLO
endif

if WITH_ROOT
ADDFFLAGS += -DWITH_ROOT
endif

if WITH_LOOPTOOLS
ADDSOURCES += LTools.F
ADDLDADD += $(LOOPTOOLS_LDFLAGS) $(LOOPTOOLS_LIBS) 
ADDFFLAGS += -DWITH_LT 

if NEW_LT
ADDFFLAGS += -DLT_NEW
endif
if OLD_LT
ADDFFLAGS += -DLT_OLD
endif
endif

if WITH_FEYNHIGGS 
ADDLDADD += $(FEYNHIGGS_LDFLAGS) $(FEYNHIGGS_LIBS) 
ADDFFLAGS += -DWITH_FH
endif

bin_PROGRAMS = vbfnlo 
vbfnlo_SOURCES = vbfnlo_main.F90 $(ADDSOURCES)
vbfnlo_LDADD = $(top_builddir)/lib/libVBFNLO.la \
               $(top_builddir)/loops/libVBFNLOLoops.la \
               $(top_builddir)/utilities/libVBFNLOUtilities.la \
               $(ADDLDADD)
vbfnlo_FFLAGS = $(AM_FFLAGS) $(ADDFFLAGS)

bin_PROGRAMS += vbfnlo_olpread 
vbfnlo_olpread_SOURCES = vbfnlo_olpread.cxx
vbfnlo_olpread_LDADD = $(top_builddir)/utilities/libVBFNLOUtilities.la \
                       $(top_builddir)/amplitudes/libVBFNLOAmplitudes.la \
                       $(top_builddir)/loops/libVBFNLOLoops.la \
                       $(top_builddir)/helas/libHELAS.la \
                       $(top_builddir)/phasespace/libVBFNLOPhasespace.la \
                       $(top_builddir)/lib/libVBFNLO.la \
                       $(FCLIBS)
vbfnlo_olpread_FFLAGS = $(AM_FFLAGS)
if WITH_LOOPS
vbfnlo_olpread_SOURCES += FFini-dummy.F
endif
if WITH_LOOPTOOLS
vbfnlo_olpread_SOURCES += LTools-dummy.F
endif

EXTRA_DIST = LTools.F.in
BUILT_SOURCES = LTools.F


if WITH_LOOPTOOLS
LTools.F: $(srcdir)/LTools.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/LTools.F.in > LTools.F
	@touch LTools.F
else
LTools.F: $(srcdir)/LTools.F.in
	@cp $(srcdir)/LTools.F.in LTools.F
	@touch LTools.F
endif

clean-local:
	rm -f $(BUILT_SOURCES)


INPUTFILES = kk_coupl_inp.dat procinfo.dat \
             vbfnlo.dat anom_HVV.dat cuts.dat kk_input.dat \
             random.dat histograms.dat \
	     susy.dat anomV.dat sps1a.slha spin2coupl.dat

dist_pkgdata_DATA = $(INPUTFILES)

printinputfiles:
	@echo $(INPUTFILES)


