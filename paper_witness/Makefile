.PHONY: all check clean
all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
check:
	python3 verify_witness.py
clean:
	latexmk -c main.tex
