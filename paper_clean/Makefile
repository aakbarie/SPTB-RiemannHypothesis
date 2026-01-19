.PHONY: all clean view

all: main.pdf

main.pdf: main.tex parts/*.tex appendices/*.tex bib/references.bib
	pdflatex main.tex
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.fls *.fdb_latexmk

view: main.pdf
	@echo "PDF compiled successfully: main.pdf"
	@echo "Page count: $$(pdfinfo main.pdf 2>/dev/null | grep Pages | awk '{print $$2}')"
