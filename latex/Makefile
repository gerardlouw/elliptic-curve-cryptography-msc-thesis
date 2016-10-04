LATEX=pdflatex
BIBTEX=bibtex
SAGE=sage

main: ecc.tex
	$(LATEX) -shell-escape ecc.tex
	$(SAGE) ecc.sagetex.sage
	$(BIBTEX) ecc
	$(LATEX) -shell-escape ecc.tex
	$(LATEX) -shell-escape ecc.tex

clean: ecc.tex
	rm -f ecc.aux ecc.lof ecc.log ecc.lot ecc.out ecc.pdf ecc.sagetex.* ecc.toc ecc.bbl ecc.blg ecc.loa commands/*.aux
	rm -rf _minted-ecc sage-plots-for-ecc.tex
