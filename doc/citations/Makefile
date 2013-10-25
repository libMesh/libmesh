bib_files := $(wildcard *.bib)

libmesh_refs: libmesh_refs.tex $(bib_files)
	pdflatex libmesh_refs.tex
	bibtex libmesh_refs
	pdflatex libmesh_refs.tex
	pdflatex libmesh_refs.tex

# Sorted in reverse-date order.  make with e.g.
# make libmesh_html SOURCE=nine.bib
# Use --nobibsource to suppress the file with raw bibtex.
# The single dash in -nofooter is apparently required?
# -d sort by date
# -r reverse the sort
# --single produce a single page with both bibtex input and output
libmesh_html: $(SOURCE)
	TMPDIR=. bibtex2html -d -r --nodoc --single -nofooter $<

clean:
	rm -rf *.aux *.dvi *.log *~ *.ps *.bbl *.blg

.PHONY : clean
