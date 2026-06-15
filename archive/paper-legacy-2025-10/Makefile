# Build automation for SPTB manuscript
TARGET = main
LATEX = latexmk
LATEXFLAGS = -pdf -interaction=nonstopmode -shell-escape

all:
	$(LATEX) $(LATEXFLAGS) $(TARGET).tex

clean:
	$(LATEX) -c $(TARGET).tex
	rm -f *.bbl *.blg *.run.xml *.out *.synctex.gz

distclean: clean
	rm -f $(TARGET).pdf
