pkglib_LTLIBRARIES = libVBFNLOPDFsets.la
libVBFNLOPDFsets_la_SOURCES = \
  Ctq6Pdf.F \
  CT10Pdf.F \
  pftopdg.F \
  mrst2004qed.F \
  mstw2008.F
PDFTABLES = \
  cteq6l1.tbl \
  ct10.00.pds \
  qed6-10gridn.dat \
  qed6-10gridp.dat \
  mstw2008lo.00.dat \
  mstw2008nlo.00.dat
pdfdatadir = $(pkgdatadir)/pdfsets
dist_pdfdata_DATA = $(PDFTABLES)


libVBFNLOPDFsets_la_FFLAGS = $(AM_FFLAGS)

$(top_builddir)/include/%.mod: libVBFNLOPDFsets_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST = Makefile.dep
