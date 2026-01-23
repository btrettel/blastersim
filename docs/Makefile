# TODO: arXiv: <https://news.ycombinator.com/item?id=29352868>
# Speed for fast feedback loop:
# - TikZ externalization: <https://tikz.dev/library-external#sec-52.4>, <https://tex.stackexchange.com/a/482560/9945>, <https://tikz.dev/base-external>
# - Ability to compile individual sections or files

# LuaLaTeX: `make TEX=lualatex`

.POSIX:

# non-POSIX
# <https://innolitics.com/articles/make-delete-on-error/>
.DELETE_ON_ERROR:
MAKEFLAGS = --warn-undefined-variables
# Later, when Make 4.4 is available, add `--shuffle` to `MAKEFLAGS`.

# <https://info.arxiv.org/help/submit_tex.html#autoproc>:
# > Naming your primary (or toplevel) file ms.tex will cause AutoTeX to process that file first in (La)TeX mode. Under PDFLaTeX modes, and in all other cases, the tex files will be processed in alphanumeric order.
KEY = blastersim

DEPS = usage.tex theory.tex verval.tex dev.tex rev.tex test_exact.tex geninput_springer.tex

.SUFFIXES: .fig .gp .tikz

BIB       = bibtex
SPELLTEX  = aspell --mode=tex --check
SPELLHTML = aspell --mode=html --check
TEX       = pdflatex
TEXFLAGS  = -halt-on-error

# Add flags like `-f` and `-l` to limit the pages if needed.
PDFTOTEXT = pdftotext

BL = \<
BR = \>

# Linux
GREP    = grep
RM      = rm
CP      = cp
DIR_SEP = /

# Windows
#GREP = findstr /n
#RM   = del
#CP   = copy
# The `^` escapes the backslash. Otherwise, `\` would be a new line in NMAKE.
#DIR_SEP = ^\

# Why `-draftmode`? This will make compilation faster as I don't want a PDF file when generated the bbl file.

$(KEY).pdf $(KEY).log: $(KEY).tex $(KEY).bbl $(DEPS)
	$(TEX) $(TEXFLAGS) -draftmode $(KEY)
	$(TEX) $(TEXFLAGS) $(KEY)

$(KEY).bbl: $(KEY).tex $(KEY).bib $(DEPS)
	$(TEX) $(TEXFLAGS) -draftmode $(KEY)
	$(BIB) $(KEY)

rev.tex: ..$(DIR_SEP)rev.tex
	$(CP) ..$(DIR_SEP)rev.tex rev.tex

test_exact.tex: ..$(DIR_SEP)test_exact.tex
	$(CP) ..$(DIR_SEP)test_exact.tex test_exact.tex

geninput_springer.tex: ..$(DIR_SEP)src$(DIR_SEP)geninput_springer.tex
	$(CP) ..$(DIR_SEP)src$(DIR_SEP)geninput_springer.tex geninput_springer.tex

# <https://math.nist.gov/~BMiller/LaTeXML/manual/commands/latexml.html>
# TODO: Fail if there are warnings. Try `--log=file`?
index.html: $(KEY).tex $(KEY).bib $(DEPS)
	latexmlc --strict --split --dest=index.html $(KEY).tex

.fig.tikz:
	fig2dev -L tikz $<

.gp.tikz:
	gnuplot $<

#$(KEY).txt: $(KEY).html
#	lynx --nolist --dump $< > $@

$(KEY).txt: $(KEY).pdf
	$(PDFTOTEXT) $(KEY).pdf $(KEY).txt

.PHONY: clean
clean:
	-$(RM) *.aux *.bbl *.blg *.css *.html *.log *.out *.toc $(KEY)_bibertool.bib

# Why not put `$(SPELLTEX) $(KEY).tex` in `$(KEY).bbl`? I get the error "Error: Stdin not a terminal." if I do that.
# `lacheck $(KEY).tex` has lots of false positives.
# TODO: Check log files (both LaTeX and BibTeX) for warnings.
# TODO: Validate HTML and CSS
# <https://github.com/dspinellis/latex-advice#latex-formatting>: > Users can define their own rules, too, for example to enforce consistency in the spelling of certain phrases. Now it is up to you to create ChkTeX rules to enforce the guidelines in this document. If you do, please share them.
# TODO: Try other BibTeX linters.
# TODO: <https://github.com/sylvainhalle/textidote>
# TODO: /home/ben/git/bmtreport/diction/
.PHONY: lint
lint: $(KEY).tex $(KEY).bib $(KEY).html $(KEY).txt $(KEY).log
	-$(GREP) ?? $(KEY).txt
	-$(GREP) "$(BL)Warning$(BR)" $(KEY).log
	$(SPELLTEX) $(KEY).tex
	$(SPELLHTML) $(KEY).html
	chktex --nowarn 2 $(KEY).tex
	lacheck $(KEY).tex
	biber --tool --validate-datamodel --dieondatamodel --quiet $(KEY).bib
	#diction --beginner --suggest $(KEY).txt
	#style --print-ari 20 $(KEY).txt
	linkchecker --check-extern --ignore-url "doi\.org/10\.1002" --ignore-url "doi\.org/10\.1061" $(KEY).html
