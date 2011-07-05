INSTALL_TARGETS = install-doc

all all-local: tutorial.html

install install-data-local: $(INSTALL_TARGETS)

tutorial.html: tutorial.xml
	$(XSLTPROC) http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl tutorial.xml >tutorial.html

EXTRA_DIST = \
	tutorial.html

MAINTAINERCLEANFILES = \
	tutorial.html
