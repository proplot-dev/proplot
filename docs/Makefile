# Minimal makefile for Sphinx documentation
# You can set these variables from the command line.
# WARNING: Sphinx 2.0 fails with numpydoc style strings so far!
# See issue: https://github.com/readthedocs/sphinx_rtd_theme/issues/750
# BUILDDIR    = ../../proplot-doc # from gh-pages workflow
SPHINXOPTS  =
SPHINXBUILD = sphinx-build
SPHINXPROJ  = ProPlot
SOURCEDIR   = .
BUILDDIR    = _build

.PHONY: help clean Makefile

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

# Make make clean ignore .git folder
# The -f doesn't raise error when files/folders not found
clean:
	rm -rf api/*
	rm -rf "$(BUILDDIR)"/html/*
	rm -rf "$(BUILDDIR)"/doctrees/*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option. See: https://github.com/sphinx-doc/sphinx/issues/6603
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)
