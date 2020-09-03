bib_files := $(filter-out templates.bib, $(wildcard *.bib))
pdf_bib_files := $(patsubst %.bib, %.pdf, $(bib_files))
LATEXMK_OPTIONS = -logfilewarninglist -g -pdf

# This target depends on all other PDFs existing.  It just combines
# all of them into a single document.
all_pdfs: all_pdfs.tex $(pdf_bib_files)
	latexmk $(LATEXMK_OPTIONS) all_pdfs.tex

# Create a single bibliography with all references.  This can be useful
# for getting a total count of all citations.
libmesh_refs: libmesh_refs.tex $(bib_files)
	latexmk $(LATEXMK_OPTIONS) libmesh_refs.tex

# Sorted in reverse-date order.  make with e.g.
# make libmesh_html SOURCE=nine.bib
# Use --nobibsource to suppress the file with raw bibtex.
# The single dash in -nofooter is apparently required?
# -d sort by date
# -r reverse the sort
# --single produce a single page with both bibtex input and output
libmesh_html: $(SOURCE)
	TMPDIR=. bibtex2html -d -r --nodoc --single -nofooter $<

# Rule for making individual year pdfs (e.g. twenty.pdf)
%.pdf: %.tex %.bib
	latexmk $(LATEXMK_OPTIONS) $<

clean:
	latexmk -C -bibtex
	rm -f *~

.PHONY : clean
