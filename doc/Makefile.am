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
	$(XSLTPROC) boostbook-xsl/docbook.xsl reference-boost.xml >reference.xml

reference-boost.xml: reference-doxy.xml boostbook-xsl
	$(XSLTPROC) boostbook-xsl/doxygen/doxygen2boostbook.xsl reference-doxy.xml >reference-boost.xml

reference-doxy.xml: xml boostbook-xsl
	$(XSLTPROC) boostbook-xsl/doxygen/collect.xsl xml/index.xml >reference-doxy.xml

xml:
	$(DOXYGEN)

# Temporary solution until we get this thing online somewhere
# boostbook-xsl: boost_1_47_0.tar.bz2
# 	tar jxvf boost_1_47_0.tar.bz2 boost_1_47_0/tools/boostbook/xsl/
# 	mv boost_1_47_0/tools/boostbook/xsl/ boostbook-xsl
# 	rm -rf boost_1_47_0

# boost_1_47_0.tar.bz2:
# 	wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2

EXTRA_DIST = \
	tutorial.xml	\
	tutorial.html

MAINTAINERCLEANFILES = \
	tutorial.html		\
	reference.html		\
	reference.xml		\
	reference-boost.xml	\
	reference-doxy.xml
