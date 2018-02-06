SUBDIRS =

pkglib_LTLIBRARIES = libVBFNLOLoops.la

libVBFNLOLoops_la_FFLAGS = $(AM_FFLAGS)
libVBFNLOLoops_la_LIBADD = 
libVBFNLOLoops_la_SOURCES = \
aux_loops.F

if WITH_LOOPS

if WITH_QCDVV
SUBDIRS += Diagrams/HexBox/aa/Dire
SUBDIRS += Diagrams/HexBox/aq/Dire
SUBDIRS += Diagrams/HexBox/aa/Cross
SUBDIRS += Diagrams/HexBox/aq/Cross
libVBFNLOLoops_la_LIBADD += Diagrams/HexBox/aa/Dire/libDiagramsAADire.la
libVBFNLOLoops_la_LIBADD += Diagrams/HexBox/aq/Dire/libDiagramsAQDire.la
libVBFNLOLoops_la_LIBADD += Diagrams/HexBox/aa/Cross/libDiagramsAACross.la
libVBFNLOLoops_la_LIBADD += Diagrams/HexBox/aq/Cross/libDiagramsAQCross.la
libVBFNLOLoops_la_FFLAGS += -DWITH_QCDVV
endif

if WITH_QCD_JJ_PLUSX
SUBDIRS += Diagrams/PenBox/
libVBFNLOLoops_la_LIBADD += Diagrams/PenBox/libDiagramsPenBox.la
endif

SUBDIRS += Diagrams/Penline
SUBDIRS += Diagrams/Boxline
libVBFNLOLoops_la_LIBADD += Diagrams/Penline/libDiagramsPenL.la
libVBFNLOLoops_la_LIBADD += Diagrams/Boxline/libDiagramsBoxL.la

SUBDIRS += TenRed/E
SUBDIRS += TenRed/D
SUBDIRS += TenRed/C
SUBDIRS += TenRed/B
libVBFNLOLoops_la_LIBADD += TenRed/E/libTIE.la
libVBFNLOLoops_la_LIBADD += TenRed/D/libTID.la
libVBFNLOLoops_la_LIBADD += TenRed/C/libTIC.la
libVBFNLOLoops_la_LIBADD += TenRed/B/libTIB.la

SUBDIRS += Integrals
libVBFNLOLoops_la_LIBADD += Integrals/libIntegrals.la

SUBDIRS += ggf_loops
libVBFNLOLoops_la_LIBADD +=  ggf_loops/libggf_loops.la

SUBDIRS += Utilities
libVBFNLOLoops_la_LIBADD +=  Utilities/libUtilities.la


if WITH_HEXAGONS
SUBDIRS += Diagrams/HexLine
libVBFNLOLoops_la_LIBADD += Diagrams/HexLine/libDiagramsHexL.la

SUBDIRS += TenRed/F
libVBFNLOLoops_la_LIBADD += TenRed/F/libTIF.la
endif

if WITH_FERMIONLOOPS
SUBDIRS += Diagrams/FermionLoops
libVBFNLOLoops_la_LIBADD += Diagrams/FermionLoops/libDiagramsFLoops.la
endif


endif

$(top_builddir)/include/%.mod: libVBFNLOLoops_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
