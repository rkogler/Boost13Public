LATEX    = latex
BIBTEX   = bibtex
DVIPS    = dvips

BASENAME = boost2013_report

default: pdflatex

nobibtex:
	latex  ${BASENAME}
	latex  ${BASENAME}
	latex  ${BASENAME}
	latex  ${BASENAME}
	dvipdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress ${BASENAME}


testlatex:
	latex  ${BASENAME}
	latex  ${BASENAME}
	bibtex ${BASENAME}
	latex  ${BASENAME}
	latex  ${BASENAME}
	dvipdf -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress ${BASENAME}

pdflatex:
	make clean
	pdflatex  ${BASENAME}
	pdflatex  ${BASENAME}
	bibtex    ${BASENAME}
	pdflatex  ${BASENAME}
	pdflatex  ${BASENAME}

pdflatexnobibtex:
	pdflatex  ${BASENAME}
	pdflatex  ${BASENAME}
	pdflatex  ${BASENAME}

#
# standard Latex targets
#

%.dvi:	%.tex 
	$(LATEX) $<

%.bbl:	%.tex *.bib
	$(LATEX) $*
	$(BIBTEX) $*

%.ps:	%.dvi
	$(DVIPS) $< -o $@

%.pdf:	%.tex
	$(PDFLATEX) $<

.PHONY: clean

clean:
	rm -f *.aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
	      *.inx *.dvi *.toc *.out *~ ~* spellTmp 

