# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

if BUILD_DOCS
DOC_TARGETS = tutorial.html
REFDOC_TARGETS = html/index.html
endif

all all-local: $(DOC_TARGETS) $(REFDOC_TARGETS)

docdir = $(datadir)/doc/lela-$(VERSION)
doc_DATA = $(DOC_TARGETS) caution.png note.png docbook.css

refdocdir = $(datadir)/doc/lela-$(VERSION)/reference
refdoc_DATA = $(REFDOC_TARGETS)

tutorial.html: $(srcdir)/config.xsl $(srcdir)/tutorial.xml
	$(XMLTO) xhtml -m $(srcdir)/config.xsl $(srcdir)/tutorial.xml
	sed -e s/index.html/tutorial.html/g index.html >tutorial.html
	rm -f index.html

html/index.html:
	$(DOXYGEN)

install-data-local:
	mkdir -p $(DESTDIR)$(refdocdir)
	if test -x $(builddir)/html; then \
		$(install_sh) $(builddir)/html/* $(DESTDIR)$(refdocdir); \
	else \
		$(install_sh) $(srcdir)/html/* $(DESTDIR)$(refdocdir); \
	fi

uninstall-hook:
	rm -rf $(DESTDIR)$(refdocdir)

maintainer-clean-local:
	rm -rf $(builddir)/html

dist-hook:
	mkdir -p $(distdir)/html
	if test -x $(builddir)/html; then \
		cp $(builddir)/html/* $(distdir)/html; \
	else \
		cp $(srcdir)/html/* $(distdir)/html; \
	fi

EXTRA_DIST = \
	tutorial.xml	\
	tutorial.html	\
	caution.png	\
	note.png	\
	config.xsl	\
	docbook.css

MAINTAINERCLEANFILES = \
	tutorial.html
