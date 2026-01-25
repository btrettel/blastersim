TEX_KEY = blastersim

TEX_DEPS = docs$(DIR_SEP)usage.tex \
docs$(DIR_SEP)theory.tex \
docs$(DIR_SEP)verval.tex \
docs$(DIR_SEP)dev.tex \
docs$(DIR_SEP)rev.tex \
docs$(DIR_SEP)test_exact.tex \
docs$(DIR_SEP)geninput_springer.tex

CLEAN_TEX = docs$(DIR_SEP)*.aux \
docs$(DIR_SEP)*.bbl \
docs$(DIR_SEP)*.blg \
docs$(DIR_SEP)*.css \
docs$(DIR_SEP)*.html \
docs$(DIR_SEP)*.log \
docs$(DIR_SEP)*.out \
docs$(DIR_SEP)*.toc \
docs$(DIR_SEP)$(TEX_KEY)_bibertool.bib

BIB        = bibtex
SPELL_TEX  = aspell --mode=tex --check
SPELL_HTML = aspell --mode=html --check
TEX        = pdflatex
TEX_FLAGS  = -halt-on-error

# Add flags like `-f` and `-l` to limit the pages if needed.
PDFTOTEXT = pdftotext

BL = \<
BR = \>

# Why `-draftmode`? This will make compilation faster as I don't want a PDF file when generated the bbl file.

docs$(DIR_SEP)$(TEX_KEY).pdf docs$(DIR_SEP)$(TEX_KEY).log: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bbl $(TEX_DEPS)
	cd docs && $(TEX) $(TEX_FLAGS) -draftmode $(TEX_KEY)
	cd docs && $(TEX) $(TEX_FLAGS) $(TEX_KEY)
	-$(GREP) "$(BL)Warning$(BR)" docs$(DIR_SEP)$(TEX_KEY).log

docs$(DIR_SEP)$(TEX_KEY).bbl: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bib $(TEX_DEPS)
	$(LOOP_START) docs$(DIR_SEP)$(TEX_KEY).tex $(TEX_DEPS) $(LOOP_MIDDLE) $(SPELL_TEX) $(LOOP_END)
	cd docs && $(TEX) $(TEX_FLAGS) -draftmode $(TEX_KEY)
	cd docs && $(BIB) $(TEX_KEY)

# Why move these files here? I don't know how to go up a directory in a cross-platform way in LaTeX. So putting everything in this directory is the easiest approach.
docs$(DIR_SEP)rev.tex: rev.tex
	$(CP) rev.tex docs$(DIR_SEP)rev.tex

docs$(DIR_SEP)test_exact.tex: test_exact.tex
	$(CP) test_exact.tex docs$(DIR_SEP)test_exact.tex

docs$(DIR_SEP)geninput_springer.tex: src$(DIR_SEP)geninput_springer.tex
	$(CP) src$(DIR_SEP)geninput_springer.tex docs$(DIR_SEP)geninput_springer.tex

# <https://math.nist.gov/~BMiller/LaTeXML/manual/commands/latexml.html>
docs$(DIR_SEP)index.html: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bib $(TEX_DEPS)
	cd docs && latexmlc --strict --split --splitnaming=label --dest=index.html $(TEX_KEY).tex
	$(LOOP_START) docs$(DIR_SEP)*.html $(LOOP_MIDDLE) $(SPELL_HTML) $(LOOP_END)
	-$(GREP) "$(BL)Warning$(BR)" docs$(DIR_SEP)$(TEX_KEY).latexml.log

docs$(DIR_SEP)$(TEX_KEY).txt: docs$(DIR_SEP)$(TEX_KEY).pdf
	cd docs && $(PDFTOTEXT) $(TEX_KEY).pdf $(TEX_KEY).txt

.PHONY: clean_tex
clean_tex:
	$(RM) $(CLEAN_TEX)
