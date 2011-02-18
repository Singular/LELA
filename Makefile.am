# Copyright (c) 2010 the LinBox group
# This file is part of LinBox
# see COPYING for licence

SUBDIRS=linbox macros tests doc

docs:doc/linbox-html/index.html

docs_dev:doc/linbox-dev-html/index.html

doc/linbox-html/index.html:
	(cd doc; make docs)

doc/linbox-dev-html/index.html:
	(cd doc; make docs_dev)

bin_SCRIPTS=linbox-config

VERSION=1.1.7
