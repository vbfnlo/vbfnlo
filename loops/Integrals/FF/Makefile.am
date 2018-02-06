noinst_LTLIBRARIES = libFF.la

if WITH_QUAD
noinst_LTLIBRARIES += libFFquad.la
endif

#Note: this list differs from the one in the Makefile because of a duplicate entry and one with an undefined function.
nodist_libFF_la_SOURCES = ff/ffinit.f
dist_libFF_la_SOURCES =  \
		ff/ff2dl2.f \
		ff/ffabcd.f \
		ff/ffca0.f \
		ff/ffcb0.f \
		ff/ffcb1.f \
		ff/ffcb2p.f \
		ff/ffcdb0.f \
		ff/ffcc0.f \
		ff/ffcc0p.f \
		ff/ffcc1.f \
		ff/ffcel2.f \
		ff/ffcel3.f \
		ff/ffcel4.f \
		ff/ffcel5.f \
		ff/ffceta.f \
		ff/ffcli2.f \
		ff/ffcrr.f \
		ff/ffcxr.f \
		ff/ffcxs3.f \
		ff/ffcxs4.f \
		ff/ffcxyz.f \
		ff/ffdcc0.f \
		ff/ffdcxs.f \
		ff/ffdel2.f \
		ff/ffdel3.f \
		ff/ffdel4.f \
		ff/ffdel5.f \
		ff/ffdel6.f \
		ff/ffdl2i.f \
		ff/ffdl5p.f \
		ff/ffdxc0.f \
		ff/ffrcvr.f \
		ff/fftran.f \
		ff/ffxb0.f \
		ff/ffxb1.f \
		ff/ffxb2p.f \
		ff/ffxc0.f \
		ff/ffxc0i.f \
		ff/ffxc0p.f \
		ff/ffxc1.f \
		ff/ffxd0.f \
		ff/ffxd0h.f \
		ff/ffxd0i.f \
		ff/ffxd0p.f \
		ff/ffxd1.f \
		ff/ffxdb0.f \
		ff/ffxdbd.f \
		ff/ffxdi.f \
		ff/ffxdpv.f \
		ff/ffxe0.f \
		ff/ffxe1.f \
		ff/ffxf0.f \
		ff/ffxf0h.f \
		ff/ffxli2.f \
		ff/ffxxyz.f

libFF_la_FFLAGS = $(AM_FCFLAGS) -I$(builddir)/ff 

FFdir = $(pkgdatadir)/ff
dist_FF_DATA = ff/fferr.dat ff/ffwarn.dat 

BUILT_SOURCES = ff/ffinit.f
EXTRA_DIST = ff/ffinit.f.in ff/ff.h.in ff/aa.h.in \
             FF.html \
             convert.sh \
             ff/Makefile \
             ff/Makefile.old \
             ff/README \
             ff/aacbc.f \
             ff/aaccc.f \
             ff/aacinv.f \
             ff/aaxbx.f \
             ff/aaxcx.f \
             ff/aaxdx.f \
             ff/aaxex.f \
             ff/aaxinv.f \
             ff/ffcb2.f \
             ff/ffcd0.f \
             ff/ffcdbd.f \
             ff/ffini.f \
             ff/ffperm5.dat \
             ff/fftest.f \
             ff/npoin.f \
             ff/npointes.f \
             ff/npointes.out \
             ff/spence.f \
             ff/update \
             ffmanual.tex 

if WITH_QUAD
# same as above with generated files
nodist_libFFquad_la_SOURCES = $(nodist_libFF_la_SOURCES:.f=.quad.F) $(dist_libFF_la_SOURCES:.f=.quad.F) 

libFFquad_la_FFLAGS = $(AM_FCFLAGS) -I$(builddir)/ff $(QUADFCFLAGS)

BUILT_SOURCES += $(nodist_libFFquad_la_SOURCES)
endif

ff/ffinit.f: $(srcdir)/ff/ffinit.f.in ff/ff.h ff/aa.h
	@mkdir -p $(builddir)/ff
	@sed s%@FFPATH@%"$(pkgdatadir)/ff"% < $(srcdir)/ff/ffinit.f.in > ff/ffinit.f
	@touch ff/ffinit.f

ff/ff.h: $(srcdir)/ff/ff.h.in
	@mkdir -p $(builddir)/ff
	@cp $< $@
	@touch $@

ff/aa.h: $(srcdir)/ff/aa.h.in
	@mkdir -p $(builddir)/ff
	@cp $< $@
	@touch $@

if WITH_QUAD
ff/ff.quad.h: ff/ff.h
	@sed s%vbf%vbfquad%g < $< > $@
	@touch $@

ff/ffinit.quad.F: $(srcdir)/ff/ffinit.f.in ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" -e "s%@FFPATH@%"$(pkgdatadir)/ff"%" < $(srcdir)/ff/ffinit.f.in > $@
	@touch $@

ff/ffinit.quad.o: ff/ffinit.quad.F ff/ff.quad.h

ff/ff2dl2.quad.F: $(srcdir)/ff/ff2dl2.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ff2dl2.quad.o: ff/ff2dl2.quad.F ff/ff.quad.h

ff/ffabcd.quad.F: $(srcdir)/ff/ffabcd.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffabcd.quad.o: ff/ffabcd.quad.F ff/ff.quad.h

ff/ffca0.quad.F: $(srcdir)/ff/ffca0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffca0.quad.o: ff/ffca0.quad.F ff/ff.quad.h

ff/ffcb0.quad.F: $(srcdir)/ff/ffcb0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcb0.quad.o: ff/ffcb0.quad.F ff/ff.quad.h

ff/ffcb1.quad.F: $(srcdir)/ff/ffcb1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcb1.quad.o: ff/ffcb1.quad.F ff/ff.quad.h

ff/ffcb2p.quad.F: $(srcdir)/ff/ffcb2p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcb2p.quad.o: ff/ffcb2p.quad.F ff/ff.quad.h

ff/ffcdb0.quad.F: $(srcdir)/ff/ffcdb0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcdb0.quad.o: ff/ffcdb0.quad.F ff/ff.quad.h

ff/ffcc0.quad.F: $(srcdir)/ff/ffcc0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcc0.quad.o: ff/ffcc0.quad.F ff/ff.quad.h

ff/ffcc0p.quad.F: $(srcdir)/ff/ffcc0p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcc0p.quad.o: ff/ffcc0p.quad.F ff/ff.quad.h

ff/ffcc1.quad.F: $(srcdir)/ff/ffcc1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcc1.quad.o: ff/ffcc1.quad.F ff/ff.quad.h

ff/ffcel2.quad.F: $(srcdir)/ff/ffcel2.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcel2.quad.o: ff/ffcel2.quad.F ff/ff.quad.h

ff/ffcel3.quad.F: $(srcdir)/ff/ffcel3.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcel3.quad.o: ff/ffcel3.quad.F ff/ff.quad.h

ff/ffcel4.quad.F: $(srcdir)/ff/ffcel4.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcel4.quad.o: ff/ffcel4.quad.F ff/ff.quad.h

ff/ffcel5.quad.F: $(srcdir)/ff/ffcel5.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcel5.quad.o: ff/ffcel5.quad.F ff/ff.quad.h

ff/ffceta.quad.F: $(srcdir)/ff/ffceta.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffceta.quad.o: ff/ffceta.quad.F ff/ff.quad.h

ff/ffcli2.quad.F: $(srcdir)/ff/ffcli2.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcli2.quad.o: ff/ffcli2.quad.F ff/ff.quad.h

ff/ffcrr.quad.F: $(srcdir)/ff/ffcrr.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcrr.quad.o: ff/ffcrr.quad.F ff/ff.quad.h

ff/ffcxr.quad.F: $(srcdir)/ff/ffcxr.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcxr.quad.o: ff/ffcxr.quad.F ff/ff.quad.h

ff/ffcxs3.quad.F: $(srcdir)/ff/ffcxs3.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcxs3.quad.o: ff/ffcxs3.quad.F ff/ff.quad.h

ff/ffcxs4.quad.F: $(srcdir)/ff/ffcxs4.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcxs4.quad.o: ff/ffcxs4.quad.F ff/ff.quad.h

ff/ffcxyz.quad.F: $(srcdir)/ff/ffcxyz.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffcxyz.quad.o: ff/ffcxyz.quad.F ff/ff.quad.h

ff/ffdcc0.quad.F: $(srcdir)/ff/ffdcc0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdcc0.quad.o: ff/ffdcc0.quad.F ff/ff.quad.h

ff/ffdcxs.quad.F: $(srcdir)/ff/ffdcxs.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdcxs.quad.o: ff/ffdcxs.quad.F ff/ff.quad.h

ff/ffdel2.quad.F: $(srcdir)/ff/ffdel2.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdel2.quad.o: ff/ffdel2.quad.F ff/ff.quad.h

ff/ffdel3.quad.F: $(srcdir)/ff/ffdel3.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdel3.quad.o: ff/ffdel3.quad.F ff/ff.quad.h

ff/ffdel4.quad.F: $(srcdir)/ff/ffdel4.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdel4.quad.o: ff/ffdel4.quad.F ff/ff.quad.h

ff/ffdel5.quad.F: $(srcdir)/ff/ffdel5.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdel5.quad.o: ff/ffdel5.quad.F ff/ff.quad.h

ff/ffdel6.quad.F: $(srcdir)/ff/ffdel6.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdel6.quad.o: ff/ffdel6.quad.F ff/ff.quad.h

ff/ffdl2i.quad.F: $(srcdir)/ff/ffdl2i.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdl2i.quad.o: ff/ffdl2i.quad.F ff/ff.quad.h

ff/ffdl5p.quad.F: $(srcdir)/ff/ffdl5p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdl5p.quad.o: ff/ffdl5p.quad.F ff/ff.quad.h

ff/ffdxc0.quad.F: $(srcdir)/ff/ffdxc0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffdxc0.quad.o: ff/ffdxc0.quad.F ff/ff.quad.h

ff/ffrcvr.quad.F: $(srcdir)/ff/ffrcvr.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffrcvr.quad.o: ff/ffrcvr.quad.F ff/ff.quad.h

ff/fftran.quad.F: $(srcdir)/ff/fftran.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/fftran.quad.o: ff/fftran.quad.F ff/ff.quad.h

ff/ffxb0.quad.F: $(srcdir)/ff/ffxb0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxb0.quad.o: ff/ffxb0.quad.F ff/ff.quad.h

ff/ffxb1.quad.F: $(srcdir)/ff/ffxb1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxb1.quad.o: ff/ffxb1.quad.F ff/ff.quad.h

ff/ffxb2p.quad.F: $(srcdir)/ff/ffxb2p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxb2p.quad.o: ff/ffxb2p.quad.F ff/ff.quad.h

ff/ffxc0.quad.F: $(srcdir)/ff/ffxc0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxc0.quad.o: ff/ffxc0.quad.F ff/ff.quad.h

ff/ffxc0i.quad.F: $(srcdir)/ff/ffxc0i.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxc0i.quad.o: ff/ffxc0i.quad.F ff/ff.quad.h

ff/ffxc0p.quad.F: $(srcdir)/ff/ffxc0p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxc0p.quad.o: ff/ffxc0p.quad.F ff/ff.quad.h

ff/ffxc1.quad.F: $(srcdir)/ff/ffxc1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxc1.quad.o: ff/ffxc1.quad.F ff/ff.quad.h

ff/ffxd0.quad.F: $(srcdir)/ff/ffxd0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxd0.quad.o: ff/ffxd0.quad.F ff/ff.quad.h

ff/ffxd0h.quad.F: $(srcdir)/ff/ffxd0h.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxd0h.quad.o: ff/ffxd0h.quad.F ff/ff.quad.h

ff/ffxd0i.quad.F: $(srcdir)/ff/ffxd0i.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxd0i.quad.o: ff/ffxd0i.quad.F ff/ff.quad.h

ff/ffxd0p.quad.F: $(srcdir)/ff/ffxd0p.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxd0p.quad.o: ff/ffxd0p.quad.F ff/ff.quad.h

ff/ffxd1.quad.F: $(srcdir)/ff/ffxd1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxd1.quad.o: ff/ffxd1.quad.F ff/ff.quad.h

ff/ffxdb0.quad.F: $(srcdir)/ff/ffxdb0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxdb0.quad.o: ff/ffxdb0.quad.F ff/ff.quad.h

ff/ffxdbd.quad.F: $(srcdir)/ff/ffxdbd.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxdbd.quad.o: ff/ffxdbd.quad.F ff/ff.quad.h

ff/ffxdi.quad.F: $(srcdir)/ff/ffxdi.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxdi.quad.o: ff/ffxdi.quad.F ff/ff.quad.h

ff/ffxdpv.quad.F: $(srcdir)/ff/ffxdpv.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxdpv.quad.o: ff/ffxdpv.quad.F ff/ff.quad.h

ff/ffxe0.quad.F: $(srcdir)/ff/ffxe0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxe0.quad.o: ff/ffxe0.quad.F ff/ff.quad.h

ff/ffxe1.quad.F: $(srcdir)/ff/ffxe1.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxe1.quad.o: ff/ffxe1.quad.F ff/ff.quad.h

ff/ffxf0.quad.F: $(srcdir)/ff/ffxf0.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxf0.quad.o: ff/ffxf0.quad.F ff/ff.quad.h

ff/ffxf0h.quad.F: $(srcdir)/ff/ffxf0h.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxf0h.quad.o: ff/ffxf0h.quad.F ff/ff.quad.h

ff/ffxli2.quad.F: $(srcdir)/ff/ffxli2.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxli2.quad.o: ff/ffxli2.quad.F ff/ff.quad.h

ff/ffxxyz.quad.F: $(srcdir)/ff/ffxxyz.f ff/ff.quad.h
	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < $< > $@
	@touch $@

ff/ffxxyz.quad.o: ff/ffxxyz.quad.F ff/ff.quad.h

endif

clean-local:
	rm -f $(BUILT_SOURCES) ff/ff.quad.h ff/aa.h ff/ff.h

#
# This script generates the entries for the .quad.F files above
# ff.h and ffinit.f need special treatment
#
# #! /bin/sh
# LIST="ff/ff2dl2.f ff/ffabcd.f ff/ffca0.f ff/ffcb0.f ff/ffcb1.f ff/ffcb2p.f ff/ffcdb0.f ff/ffcc0.f ff/ffcc0p.f ff/ffcc1.f ff/ffcel2.f ff/ffcel3.f ff/ffcel4.f ff/ffcel5.f ff/ffceta.f ff/ffcli2.f ff/ffcrr.f ff/ffcxr.f ff/ffcxs3.f ff/ffcxs4.f ff/ffcxyz.f ff/ffdcc0.f ff/ffdcxs.f ff/ffdel2.f ff/ffdel3.f ff/ffdel4.f ff/ffdel5.f ff/ffdel6.f ff/ffdl2i.f ff/ffdl5p.f ff/ffdxc0.f ff/ffrcvr.f ff/fftran.f ff/ffxb0.f ff/ffxb1.f ff/ffxb2p.f ff/ffxc0.f ff/ffxc0i.f ff/ffxc0p.f ff/ffxc1.f ff/ffxd0.f ff/ffxd0h.f ff/ffxd0i.f ff/ffxd0p.f ff/ffxd1.f ff/ffxdb0.f ff/ffxdbd.f ff/ffxdi.f ff/ffxdpv.f ff/ffxe0.f ff/ffxe1.f ff/ffxf0.f ff/ffxf0h.f ff/ffxli2.f ff/ffxxyz.f"
# 
# for i in $LIST; do
# quadin=`echo $i | sed -e "s/\.f/.quad.F/"`
# quadout=`echo $i | sed -e "s/\.f/.quad.o/"`
# cat <<EOF
# $quadin: \$(srcdir)/$i ff/ff.quad.h
# 	@sed -e "s%vbf%vbfquad%g" -e "s%ff\.h%ff.quad.h%g" < \$< > \$@
# 	@touch \$@
# 
# $quadout: $quadin ff/ff.quad.h
# 
# EOF
# done

$(top_builddir)/include/%.mod: libFF_la-%.lo
	$(MAKE) $^

include $(top_srcdir)/Makefile.dependencylist
EXTRA_DIST += Makefile.dep
