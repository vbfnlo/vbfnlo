pkglib_LTLIBRARIES = libVBFNLOUtilities.la
libVBFNLOUtilities_la_SOURCES = \
process.inc \
coupl_haddecay.inc \
mssm.inc \
hist_real.inc \
coupl.inc \
tau_hel.inc \
BLHAhelper.inc \
scales.inc \
koppln.inc \
cuts.inc \
ckm.inc \
mgtest.inc \
spin2coupl.inc \
lha.inc \
Apdf.inc \
hist.inc \
an_couplings.inc \
global.inc \
kkcoupl.inc \
abbr_NC.h \
SLHADefs.h \
roothist.h \
abbr_CC.h \
PDG.h \
BLHAinterface.h \
integrals_NC.h \
SLHAReadBlocks.h \
integrals_CC.h \
SLHA.h \
ff.h \
BLHAclass.h \
SLHAIntern.h \
globalvars.F90 \
instabilities.F \
readinput.F90 \
cmd_args.F90 \
statistics.F\
brakets.F \
diptest.F \
dotproducts.F \
kinematics_aux.F \
hist_stor.F90 \
monaco_rng_mz.F90 \
monaco_rng_sob.F90 \
random/xoroshiro128plus.c \
random/splitmix64.c \
monaco.F90 \
vbfnlo_mpi.F90 \
stopwatch.F90 \
cuts.F \
choose_ps.F \
jetdefinition.F \
leshouches.F \
ptilde.F \
reweighting.F \
histograms.F \
random.F90 \
scales.F90 \
globals.F \
HtoWW.F \
anomHiggs.F \
anomV.F \
finitecoll.F \
CataniSeymour.F\
FCtools.F \
writedata.c \
ckm.F \
BLHAhelper.F90 \
BLHAmommapping.F90 \
BLHAinterface.cpp \
BLHAclass.cpp \
stdoutredirect.F90 \
CPS_HIGGSTO.F90 \
ggf_init.F \
process.F \
koppln.F \
coupl_haddecay.F \
parameters.F \
slha_read.F \
pdfsets.F90 \
vbfnloversion.F90 \
getPaths.F90



libVBFNLOUtilities_la_FFLAGS = $(AM_FFLAGS)
libVBFNLOUtilities_la_CPPFLAGS = $(AM_CPPFLAGS) 
libVBFNLOUtilities_la_CFLAGS = $(AM_CFLAGS) -std=c99

libVBFNLOUtilities_la_LIBADD= $(top_builddir)/PDFsets/libVBFNLOPDFsets.la

libVBFNLOUtilities_la_FFLAGS += -DINPUTPATH=\"'$(pkgdatadir)'\" -DPDFSETSPATH=\"'$(pkgdatadir)/pdfsets'\"

if WITH_GGF
libVBFNLOUtilities_la_FFLAGS += -DWITH_GGF
endif

if WITH_DIBOSON
libVBFNLOUtilities_la_FFLAGS += -DWITH_DIBOSON
endif

if WITH_DIBOSONJET
libVBFNLOUtilities_la_FFLAGS += -DWITH_DIBOSONJET
endif

if WITH_QCDVV
libVBFNLOUtilities_la_FFLAGS += -DWITH_QCDVV
endif

if WITH_TRIBOSON
libVBFNLOUtilities_la_FFLAGS += -DWITH_TRIBOSON
endif

if WITH_TRIBOSONJET
libVBFNLOUtilities_la_FFLAGS += -DWITH_TRIBOSONJET
endif

if WITH_QCDV
libVBFNLOUtilities_la_FFLAGS += -DWITH_QCDV
endif

if WITH_VBF
libVBFNLOUtilities_la_FFLAGS += -DWITH_VBF
endif

if WITH_HJJJ
libVBFNLOUtilities_la_FFLAGS += -DWITH_HJJJ
endif

if WITH_KK
libVBFNLOUtilities_la_FFLAGS += -DWITH_KK
endif

libVBFNLOUtilities_la_SOURCES += \
spin2width.F \
spin2coupl.F

if WITH_NLO
libVBFNLOUtilities_la_FFLAGS += -DWITH_NLO
endif

if WITH_MPI
libVBFNLOUtilities_la_FFLAGS += -DWITH_MPI
endif

if WITH_LHAPDF
libVBFNLOUtilities_la_FFLAGS += -DWITH_LHA
libVBFNLOUtilities_la_LIBADD += $(LHAPDF_LDFLAGS) $(LHAPDF_LIBS)
endif

if WITH_LHAPDF_6
libVBFNLOUtilities_la_FFLAGS += -DWITH_LHAPDF6
endif

if WITH_FEYNHIGGS
libVBFNLOUtilities_la_SOURCES += FHiggs.F
libVBFNLOUtilities_la_FFLAGS += -DWITH_FH -I $(FEYNHIGGS_DIR)/include

libVBFNLOUtilities_la_LIBADD += $(FEYNHIGGS_LDFLAGS) $(FEYNHIGGS_LIBS)

libVBFNLOUtilities_la_FFLAGS += -DFHVERSION=$(FHVERSION)

endif

if WITH_LOOPTOOLS
libVBFNLOUtilities_la_SOURCES += \
ew_virtual.F \
NC_boxes.F \
CC_boxes.F \
CC_qqV.F \
CC_pentagons.F \
NC_qqV.F \
NC_pentagons.F \
ew_formfactors.F \
ew_qqV_eff.F \
mat_elCC+.F \
renormalisation.F \
mat_elNC+.F

libVBFNLOUtilities_la_FFLAGS += -DWITH_LT -I $(LOOPTOOLS_DIR)/include

libVBFNLOUtilities_la_LIBADD += $(LOOPTOOLS_LDFLAGS) $(LOOPTOOLS_LIBS)
endif

if WITH_ROOT
libVBFNLOUtilities_la_FFLAGS += -DWITH_ROOT
libVBFNLOUtilities_la_CPPFLAGS += $(ROOTINCLUDE)
libVBFNLOUtilities_la_LIBADD += $(ROOTLDFLAGS) $(ROOTLIBS)
libVBFNLOUtilities_la_SOURCES += roothist.cpp rootuserhists.cpp
endif

if WITH_HEPMC
libVBFNLOUtilities_la_FFLAGS += -DWITH_HEPMC
libVBFNLOUtilities_la_CPPFLAGS += $(HEPMCINCLUDE)
libVBFNLOUtilities_la_LIBADD += $(HEPMCLDFLAGS) $(HEPMCLIBS)
libVBFNLOUtilities_la_SOURCES += hepmc.cpp 
endif

BUILT_SOURCES = vbfnloversion.inc

vbfnloversion.inc:
	@echo strgithash = \"GIT `git rev-parse --short=8 HEAD`\" > vbfnloversion.inc

.PHONY: vbfnloversion.inc

clean-local:
	rm -f $(BUILT_SOURCES) hist_stor.mod

if WITH_KK
libVBFNLOUtilities_la_FFLAGS += -DWITH_KK
libVBFNLOUtilities_la_CPPFLAGS += $(GSLINCLUDE)
libVBFNLOUtilities_la_LIBADD += $(GSLLIBS)
libVBFNLOUtilities_la_SOURCES += \
kkcoupl.cpp \
kk_coupl_inp.F  \
mhlm_coupl.F \
mhlm_functions.F
endif


if WITH_QUAD
libVBFNLOUtilities_la_FFLAGS += -DWITH_QUAD
endif


##
## for the .F90 files libtool uses FCFLAGS instead of FFLAGS
## => make the two equal, since up to now only FFLAGS has been
##    extended in the .am file...
##
libVBFNLOUtilities_la_FCFLAGS = $(libVBFNLOUtilities_la_FFLAGS)

$(top_builddir)/include/%.mod: libVBFNLOUtilities_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
