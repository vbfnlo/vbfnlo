BUILT_SOURCES = done-all-links

AUTOMAKE_OPTIONS = -Wno-portability

DIRLINKS = $(abs_top_srcdir)/amplitudes \
		   $(abs_top_srcdir)/helas \
		   $(abs_top_srcdir)/loops \
           $(abs_top_srcdir)/PDFsets \
           $(abs_top_srcdir)/phasespace \
           $(abs_top_srcdir)/utilities

CLEANFILES = done-all-links 

done-all-links: $(DIRLINKS)
	mkdir -p VBFNLO
	$(LN_S) -f $(DIRLINKS) VBFNLO
	touch done-all-links

install-data-local:
	find VBFNLO -follow \( -name '*.h' -or -name '*.inc' \) \
          -exec $(install_sh_DATA) \{\} $(DESTDIR)$(includedir)/\{\} \;

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/VBFNLO

clean-local:
	rm -rf VBFNLO
	rm -f *.mod

