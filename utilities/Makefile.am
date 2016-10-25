pkglib_LTLIBRARIES = libVBFNLOUtilities.la libVBFspecific.la libGGFspecific.la
libVBFNLOUtilities_la_SOURCES = \
instabilities.F\
statistics.F\
brakets.F \
diptest.F \
dotproducts.F \
kinematics_aux.F \
monaco8.F \
readinput.F \
getPaths.F90 \
time.F \
cuts.F \
choose_ps.F \
jetdefinition.F \
leshouches.F \
ptilde.F \
reweighting.F \
histograms.F \
histogram_storage.F90 \
random.F \
scales.F \
globals.F \
HtoWW.F \
anomHiggs.F \
anomV.F \
pdfsets.F90 \
finitecoll.F \
CataniSeymour.F\
FCtools.F \
version.F \
ckm.F \
BLHAhelper.F90 \
BLHAmommapping.F90 \
BLHAinterface.cpp \
BLHAclass.cpp \
CPS_HIGGSTO.F

libVBFspecific_la_SOURCES = \
ggf_init.F \
process.F \
koppln.F \
coupl_haddecay.F \
parameters.F \
slha_read.F 

libGGFspecific_la_SOURCES = \
ggf_init.F \
process.F \
koppln.F \
coupl_haddecay.F \
parameters.F \
slha_read.F 

libVBFNLOUtilities_la_FFLAGS = $(AM_FFLAGS)
libVBFNLOUtilities_la_CPPFLAGS = $(AM_CPPFLAGS) 

libVBFNLOUtilities_la_LIBADD= \
$(top_builddir)/PDFsets/libVBFNLOPDFsets.la

libVBFspecific_la_FFLAGS = $(AM_FFLAGS)
libVBFspecific_la_LIBADD= \
libVBFNLOUtilities.la

libGGFspecific_la_FFLAGS = $(AM_FFLAGS)
libGGFspecific_la_LIBADD= \
libVBFNLOUtilities.la


if WITH_GGF
libGGFspecific_la_FFLAGS += -DWITH_GGF
endif

if WITH_DIBOSON
libVBFspecific_la_FFLAGS += -DWITH_DIBOSON
libGGFspecific_la_FFLAGS += -DWITH_DIBOSON
endif

if WITH_DIBOSONJET
libVBFspecific_la_FFLAGS += -DWITH_DIBOSONJET
libGGFspecific_la_FFLAGS += -DWITH_DIBOSONJET
endif

if WITH_QCDVV
libVBFspecific_la_FFLAGS += -DWITH_QCDVV
endif

if WITH_TRIBOSON
libVBFspecific_la_FFLAGS += -DWITH_TRIBOSON
endif

if WITH_TRIBOSONJET
libVBFspecific_la_FFLAGS += -DWITH_TRIBOSONJET
endif

if WITH_QCDV
libVBFspecific_la_FFLAGS += -DWITH_QCDV
endif

if WITH_VBF
libVBFNLOUtilities_la_FFLAGS += -DWITH_VBF
libVBFspecific_la_FFLAGS += -DWITH_VBF
endif

if WITH_HJJJ
libVBFNLOUtilities_la_FFLAGS += -DWITH_HJJJ
libVBFspecific_la_FFLAGS += -DWITH_HJJJ
endif

if WITH_KK
libVBFspecific_la_FFLAGS += -DWITH_KK
endif

if WITH_SPIN2
libVBFspecific_la_SOURCES += \
spin2width.F \
spin2coupl.F
libVBFspecific_la_FFLAGS += -DWITH_SPIN2
endif

if WITH_NLO
libVBFNLOUtilities_la_FFLAGS += -DWITH_NLO
libVBFspecific_la_FFLAGS += -DWITH_NLO
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
libVBFNLOUtilities_la_FFLAGS += -DWITH_FH
libVBFspecific_la_FFLAGS += -DWITH_FH

if WITH_GGFLO
libGGFspecific_la_FFLAGS += -DWITH_FH
endif

libVBFNLOUtilities_la_LIBADD += $(FEYNHIGGS_LDFLAGS) $(FEYNHIGGS_LIBS)

if FH26
libVBFNLOUtilities_la_FFLAGS += -DWITH_FH26
endif
if FH27
libVBFNLOUtilities_la_FFLAGS += -DWITH_FH27
endif
if FH28
libVBFNLOUtilities_la_FFLAGS += -DWITH_FH28
endif
if FH_SHARED
libVBFNLOUtilities_la_FFLAGS += -DFH_SHRD
endif

endif

if WITH_LOOPTOOLS
libVBFNLOUtilities_la_SOURCES += \
renormalisation.F \
ew_formfactors.F \
ew_qqV_eff.F \
ew_virtual.F \
CC_boxes.F \
CC_qqV.F \
CC_pentagons.F \
mat_elCC+.F \
NC_boxes.F \
NC_qqV.F \
NC_pentagons.F \
mat_elNC+.F
libVBFNLOUtilities_la_FFLAGS += -DWITH_LT
libVBFspecific_la_FFLAGS += -DWITH_LT
if WITH_GGFLO
libGGFspecific_la_FFLAGS += -DWITH_LT
endif
libVBFNLOUtilities_la_LIBADD += $(LOOPTOOLS_LDFLAGS) $(LOOPTOOLS_LIBS)
endif

if WITH_ROOT
libVBFNLOUtilities_la_FFLAGS += -DWITH_ROOT
libVBFNLOUtilities_la_CPPFLAGS += $(ROOTINCLUDE)
libVBFNLOUtilities_la_LIBADD += $(ROOTLDFLAGS) $(ROOTLIBS)
libVBFNLOUtilities_la_SOURCES += \
rootevents.cpp \
roothist.cpp \
rootuserhists.cpp
endif

if WITH_HEPMC
libVBFNLOUtilities_la_FFLAGS += -DWITH_HEPMC
libVBFNLOUtilities_la_CPPFLAGS += $(HEPMCINCLUDE)
libVBFNLOUtilities_la_LIBADD += $(HEPMCLDFLAGS) $(HEPMCLIBS)
libVBFNLOUtilities_la_SOURCES += \
hepmc.cpp 
endif

EXTRA_DIST = getPaths.F90.in pdfsets.F90.in FHiggs.F.in renormalisation.F.in ew_qqV_eff.F.in ew_formfactors.F.in mat_elCC+.F.in CC_pentagons.F.in CC_boxes.F.in CC_qqV.F.in mat_elNC+.F.in NC_pentagons.F.in NC_boxes.F.in NC_qqV.F.in version.F.in
BUILT_SOURCES = getPaths.F90 pdfsets.F90 FHiggs.F renormalisation.F ew_qqV_eff.F ew_formfactors.F mat_elCC+.F CC_pentagons.F CC_boxes.F CC_qqV.F mat_elNC+.F NC_pentagons.F NC_boxes.F NC_qqV.F version.F
if USE_SVNVERSION
BUILT_SOURCES += version.tmp
endif


## extra dependencies for the histogram routine: the include files (which are only used in histograms.F) and the
## module file hist_store.mod from histogram_storage.F90 (which is a by-product of creating the .lo file)
## hist_stor.mod is cleaned by make clean explicitly by adding it to 'clean-local'
histograms.F: $(srcdir)/hist.inc $(srcdir)/hist_real.inc libVBFNLOUtilities_la-histogram_storage.lo
	@touch $(srcdir)/histograms.F

getPaths.F90: $(srcdir)/getPaths.F90.in
	@sed -e s%@INPUTPATH@%"$(pkgdatadir)"% -e s%@PDFSETSPATH@%"$(pkgdatadir)/pdfsets"% < $(srcdir)/getPaths.F90.in > getPaths.F90
	@touch getPaths.F90

if WITH_LHAPDF
pdfsets.F90: $(srcdir)/pdfsets.F90.in
	@sed s%@LHAPDFPATH@%"$(LHAPDF_PDFPATH)"% < $(srcdir)/pdfsets.F90.in > pdfsets.F90
	@touch pdfsets.F90
else
pdfsets.F90: $(srcdir)/pdfsets.F90.in
	@cp $(srcdir)/pdfsets.F90.in pdfsets.F90
	@touch pdfsets.F90
endif


if WITH_FEYNHIGGS
FHiggs.F: $(srcdir)/FHiggs.F.in
	@sed s%@FEYNHIGGSPATH@%"$(FEYNHIGGS_DIR)"% < $(srcdir)/FHiggs.F.in > FHiggs.F
	@touch FHiggs.F
else
FHiggs.F: $(srcdir)/FHiggs.F.in
	@cp $(srcdir)/FHiggs.F.in FHiggs.F
	@touch FHiggs.F
endif


if WITH_LOOPTOOLS
renormalisation.F: $(srcdir)/renormalisation.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/renormalisation.F.in > renormalisation.F
	@touch renormalisation.F

ew_qqV_eff.F: $(srcdir)/ew_qqV_eff.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/ew_qqV_eff.F.in > ew_qqV_eff.F
	@touch ew_qqV_eff.F

ew_formfactors.F: $(srcdir)/ew_formfactors.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/ew_formfactors.F.in > ew_formfactors.F
	@touch ew_formfactors.F

mat_elCC+.F: $(srcdir)/mat_elCC+.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/mat_elCC+.F.in > mat_elCC+.F
	@touch mat_elCC+.F

CC_pentagons.F: $(srcdir)/CC_pentagons.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/CC_pentagons.F.in > CC_pentagons.F
	@touch CC_pentagons.F

CC_boxes.F: $(srcdir)/CC_boxes.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/CC_boxes.F.in > CC_boxes.F
	@touch CC_boxes.F

CC_qqV.F: $(srcdir)/CC_qqV.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/CC_qqV.F.in > CC_qqV.F
	@touch CC_qqV.F

mat_elNC+.F: $(srcdir)/mat_elNC+.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/mat_elNC+.F.in > mat_elNC+.F
	@touch mat_elNC+.F

NC_pentagons.F: $(srcdir)/NC_pentagons.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/NC_pentagons.F.in > NC_pentagons.F
	@touch NC_pentagons.F

NC_boxes.F: $(srcdir)/NC_boxes.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/NC_boxes.F.in > NC_boxes.F
	@touch NC_boxes.F

NC_qqV.F: $(srcdir)/NC_qqV.F.in
	@sed s%@LOOPTOOLSPATH@%"$(LOOPTOOLS_DIR)"% < $(srcdir)/NC_qqV.F.in > NC_qqV.F
	@touch NC_qqV.F

else
renormalisation.F: $(srcdir)/renormalisation.F.in
	@cp $(srcdir)/renormalisation.F.in renormalisation.F
	@touch renormalisation.F

ew_qqV_eff.F: $(srcdir)/ew_qqV_eff.F.in
	@cp $(srcdir)/ew_qqV_eff.F.in ew_qqV_eff.F
	@touch ew_qqV_eff.F

ew_formfactors.F: $(srcdir)/ew_formfactors.F.in
	@cp $(srcdir)/ew_formfactors.F.in ew_formfactors.F
	@touch ew_formfactors.F

mat_elCC+.F: $(srcdir)/mat_elCC+.F.in
	@cp $(srcdir)/mat_elCC+.F.in mat_elCC+.F
	@touch mat_elCC+.F

CC_pentagons.F: $(srcdir)/CC_pentagons.F.in
	@cp $(srcdir)/CC_pentagons.F.in CC_pentagons.F
	@touch CC_pentagons.F

CC_boxes.F: $(srcdir)/CC_boxes.F.in
	@cp $(srcdir)/CC_boxes.F.in CC_boxes.F
	@touch CC_boxes.F

CC_qqV.F: $(srcdir)/CC_qqV.F.in
	@cp $(srcdir)/CC_qqV.F.in CC_qqV.F
	@touch CC_qqV.F

mat_elNC+.F: $(srcdir)/mat_elNC+.F.in
	@cp $(srcdir)/mat_elNC+.F.in mat_elNC+.F
	@touch mat_elNC+.F

NC_pentagons.F: $(srcdir)/NC_pentagons.F.in
	@cp $(srcdir)/NC_pentagons.F.in NC_pentagons.F
	@touch NC_pentagons.F

NC_boxes.F: $(srcdir)/NC_boxes.F.in
	@cp $(srcdir)/NC_boxes.F.in NC_boxes.F
	@touch NC_boxes.F

NC_qqV.F: $(srcdir)/NC_qqV.F.in
	@cp $(srcdir)/NC_qqV.F.in NC_qqV.F
	@touch NC_qqV.F
endif

if USE_SVNVERSION
.PHONY: version
version:
	@(( test -d $(top_srcdir)/.svn && svnversion -cn $(top_srcdir)) \
		|| (test -d $(top_srcdir)/.git && git --no-pager log -1 --pretty=format:%h )) \
         | sed -e "s/.*://" -e "s/\([0-9a-f]*\).*/\1/" | grep "[0-9a-f]" > version.tmp.new

version.tmp: version
	@if ! diff -q version.tmp{.new,}; then mv -f version.tmp{.new,}; fi
	@rm -f version.tmp.new

version.F: $(srcdir)/version.F.in version.tmp
	@sed s%@VBFNLOSVNVERSION@%`cat version.tmp`% $(srcdir)/version.F.in > version.F
	@touch version.F
else 
version.F: $(srcdir)/version.F.in
	@sed s%@VBFNLOSVNVERSION@%% \
          $(srcdir)/version.F.in > version.F
	@touch version.F
endif 

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

