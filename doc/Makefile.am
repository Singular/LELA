INSTALL_TARGETS = install-doc

if BUILD_DOCS
DOC_TARGETS = tutorial.html
endif

all all-local: $(DOC_TARGETS)

install install-data-local: $(INSTALL_TARGETS)

tutorial.html: tutorial.xml
	$(XSLTPROC) http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl tutorial.xml >tutorial.html

EXTRA_DIST = \
	tutorial.html

MAINTAINERCLEANFILES = \
	tutorial.html
