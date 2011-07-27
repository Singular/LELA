# Copyright 2011 Bradford Hovinen
# 
# This file is part of LELA, licensed under the GNU General Public
# License version 3. See COPYING for more information.

if BUILD_DOCS
DOC_TARGETS = tutorial.html
REFDOC_TARGETS = 
endif

all all-local: $(DOC_TARGETS) $(REFDOC_TARGETS)

docdir = $(datadir)/doc/lela-$(VERSION)
doc_DATA = $(DOC_TARGETS)

refdocdir = $(datadir)/doc/lela-$(VERSION)/reference
refdoc_DATA = $(REFDOC_TARGETS)

tutorial.html: tutorial.xml
	$(XMLTO) xhtml -m config.xsl tutorial.xml
	sed -e s/index.html/tutorial.html/g index.html >tutorial.html
	rm -f index.html

html:
	$(DOXYGEN)

EXTRA_DIST = \
	tutorial.xml	\
	tutorial.html	\
	caution.png	\
	note.png	\
	config.xsl	\
	docbook.css

MAINTAINERCLEANFILES = \
	tutorial.html
