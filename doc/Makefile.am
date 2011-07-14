if BUILD_DOCS
DOC_TARGETS = tutorial.html reference.html
else
DOC_TARGETS =
endif

all all-local: $(DOC_TARGETS)

docdir = $(datadir)/doc/lela-$(VERSION)
doc_DATA = $(DOC_TARGETS)

tutorial.html: tutorial.xml
	$(XSLTPROC) http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl tutorial.xml >tutorial.html

reference.html: reference.xml
	$(XSLTPROC) http://docbook.sourceforge.net/release/xsl/current/html/docbook.xsl reference.xml >reference.html

reference.xml: reference-boost.xml boostbook-xsl
	$(XSLTPROC) boostbook-xsl/docbook.xsl reference-boost.xml >reference.xml

reference-boost.xml: reference-doxy.xml boostbook-xsl
	$(XSLTPROC) boostbook-xsl/doxygen/doxygen2boostbook.xsl reference-doxy.xml >reference-boost.xml

reference-doxy.xml: xml boostbook-xsl
	$(XSLTPROC) boostbook-xsl/doxygen/collect.xsl xml/index.xml >reference-doxy.xml

xml:
	$(DOXYGEN)

# Temporary solution until we get this thing online somewhere
boostbook-xsl: boost_1_47_0.tar.bz2
	tar jxvf boost_1_47_0.tar.bz2 boost_1_47_0/tools/boostbook/xsl/
	mv boost_1_47_0/tools/boostbook/xsl/ boostbook-xsl
	rm -rf boost_1_47_0

boost_1_47_0.tar.bz2:
	wget http://downloads.sourceforge.net/project/boost/boost/1.47.0/boost_1_47_0.tar.bz2

EXTRA_DIST = \
	tutorial.html		\
	reference.html

MAINTAINERCLEANFILES = \
	tutorial.html		\
	reference.html		\
	reference.xml		\
	reference-boost.xml	\
	reference-doxy.xml
