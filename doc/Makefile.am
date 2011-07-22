if BUILD_DOCS
DOC_TARGETS = tutorial.html
REFDOC_TARGETS = 
endif

all all-local: $(DOC_TARGETS) $(REFDOC_TARGETS)

docdir = $(datadir)/doc/lela-$(VERSION)
doc_DATA = $(DOC_TARGETS)

refdocdir = $(datadir)/doc/lela-$(VERSION)/reference

tutorial.html: tutorial.xml
	$(XSLTPROC) http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl tutorial.xml >tutorial.html

reference: reference.xml
	mkdir reference
	$(XSLTPROC) --stringparam base.dir reference/ http://docbook.sourceforge.net/release/xsl/current/html/chunk.xsl reference.xml

reference.xml: reference-boost.xml boostbook-xsl
	$(XSLTPROC) http://www.singular.uni-kl.de/lela/boostbook-xsl/docbook.xsl reference-boost.xml >reference.xml

reference-boost.xml: reference-doxy.xml boostbook-xsl
	$(XSLTPROC) http://www.singular.uni-kl.de/lela/boostbook-xsl/doxygen/doxygen2boostbook.xsl reference-doxy.xml >reference-boost.xml

reference-doxy.xml: xml boostbook-xsl
	$(XSLTPROC) http://www.singular.uni-kl.de/lela/boostbook-xsl/doxygen/collect.xsl xml/index.xml >reference-doxy.xml

xml:
	$(DOXYGEN)

EXTRA_DIST = \
	tutorial.xml	\
	tutorial.html

MAINTAINERCLEANFILES = \
	tutorial.html		\
	reference.html		\
	reference.xml		\
	reference-boost.xml	\
	reference-doxy.xml
