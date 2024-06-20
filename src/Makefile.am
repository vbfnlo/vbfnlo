bin_PROGRAMS = vbfnlo 
vbfnlo_SOURCES = vbfnlo_main.F90
vbfnlo_LDADD = $(top_builddir)/lib/libVBFNLO.la \
               $(top_builddir)/loops/libVBFNLOLoops.la \
               $(top_builddir)/utilities/libVBFNLOUtilities.la
vbfnlo_FFLAGS = $(AM_FFLAGS)

if WITH_LOOPS
vbfnlo_SOURCES += FFini.F90
if WITH_QUAD
vbfnlo_FFLAGS += -DWITH_QUAD
endif
endif

if WITH_NLO
vbfnlo_FFLAGS += -DWITH_NLO
endif

if WITH_ROOT
vbfnlo_FFLAGS += -DWITH_ROOT
endif

if WITH_LOOPTOOLS
vbfnlo_SOURCES += LTools.F
vbfnlo_LDADD += $(LOOPTOOLS_LDFLAGS) $(LOOPTOOLS_LIBS) 
vbfnlo_FFLAGS += -DWITH_LT -I $(LOOPTOOLS_DIR)/include

if NEW_LT
vbfnlo_FFLAGS += -DLT_NEW
endif
if OLD_LT
vbfnlo_FFLAGS += -DLT_OLD
endif
endif

if WITH_FEYNHIGGS 
vbfnlo_LDADD += $(FEYNHIGGS_LDFLAGS) $(FEYNHIGGS_LIBS) 
vbfnlo_FFLAGS += -DWITH_FH
endif

clean-local:
	rm -f $(BUILT_SOURCES)


INPUTFILES = kk_coupl_inp.dat procinfo.dat \
             vbfnlo.dat anom_HVV.dat cuts.dat kk_input.dat \
             random.dat histograms.dat ggflo.dat \
             susy.dat anomV.dat sps1a.slha spin2coupl.dat \
             LHC-Matchbox.in

dist_pkgdata_DATA = $(INPUTFILES)

printinputfiles:
	@echo $(INPUTFILES)

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
