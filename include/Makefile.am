BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(top_srcdir)/amplitudes \
	   $(top_srcdir)/helas \
	   $(top_srcdir)/loops \
           $(top_srcdir)/PDFsets \
           $(top_srcdir)/phasespace \
           $(top_srcdir)/utilities

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS)
	mkdir -p VBFNLO
	$(LN_S) -f $(addprefix ../, $(DIRLINKS)) VBFNLO
	touch done-all-links

install-data-local:
	find VBFNLO -follow \( -name '*.h' -or -name '*.inc' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/VBFNLO

clean-local:
	rm -rf VBFNLO

