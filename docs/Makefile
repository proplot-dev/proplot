# Minimal makefile for Sphinx documentation
# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = ProPlot
SOURCEDIR     = .
BUILDDIR      = ../../proplot-doc
# BUILDDIR      = build

.PHONY: help clean Makefile

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Make make clean ignore .git folder
# The -f doesn't raise error when files/folders not found
clean:
	rm -f api/*
	rm -f _autosummary/*
	rm -rf "$(BUILDDIR)"/html/*
	rm -rf "$(BUILDDIR)"/doctrees/*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

