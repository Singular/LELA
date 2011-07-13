# Copyright (c) 2010 the LELA group
# This file is part of LELA
# see COPYING for licence

SUBDIRS=lela util macros tests doc

docs:doc/lela-html/index.html

docs_dev:doc/lela-dev-html/index.html

doc/lela-html/index.html:
	(cd doc; make docs)

doc/lela-dev-html/index.html:
	(cd doc; make docs_dev)

bin_SCRIPTS=lela-config

VERSION=1.1.7
