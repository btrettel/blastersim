TEX_KEY = blastersim

SPELL_DEPS = docs$(DIR_SEP)usage.tex \
docs$(DIR_SEP)theory.tex \
docs$(DIR_SEP)verval.tex \
docs$(DIR_SEP)dev.tex \
docs$(DIR_SEP)springer-figures.tex \
docs$(DIR_SEP)exact-solution-figure.tex

TEX_GEN = docs$(DIR_SEP)rev.tex \
docs$(DIR_SEP)test_exact.tex \
docs$(DIR_SEP)geninput_springer.tex \
docs$(DIR_SEP)springer-example.nml \
docs$(DIR_SEP)blastersim-out-1.txt \
docs$(DIR_SEP)blastersim-out-2.txt \
docs$(DIR_SEP)defaults.tex

TEX_DEPS = mk$(DIR_SEP)latex.mk \
$(SPELL_DEPS) \
$(TEX_GEN) \

CLEAN_TEX = $(TEX_GEN) \
docs$(DIR_SEP)*.aux \
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
	-$(GREP) "$(BL)pdfTeX warning$(BR)" docs$(DIR_SEP)$(TEX_KEY).log
	-$(GREP) "$(BL)Overfull \\\hbox$(BR)" docs$(DIR_SEP)$(TEX_KEY).log
	-$(GREP) "$(BL)Rerun$(BR)" docs$(DIR_SEP)$(TEX_KEY).blg
	-$(GREP) "$(BL)Warning$(BR)" docs$(DIR_SEP)$(TEX_KEY).blg

# ChkTeX doesn't handle `\input` macros correctly unless it's in the same directory as the TeX files. The `TeXInputs` options of ChkTeX seems to use only absolute directory, so that's not a workaround.
docs$(DIR_SEP)$(TEX_KEY).bbl: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bib $(TEX_DEPS)
	$(PYTHON) py$(DIR_SEP)tripwire.py $(ALLSRC)
	$(LOOP_START) docs$(DIR_SEP)$(TEX_KEY).tex $(SPELL_DEPS) $(LOOP_MIDDLE) $(SPELL_TEX) $(LOOP_END)
	cd docs && chktex --localrc chktexrc $(TEX_KEY).tex
	biber --tool --validate-datamodel --dieondatamodel --quiet docs$(DIR_SEP)$(TEX_KEY).bib
	cd docs && $(TEX) $(TEX_FLAGS) -draftmode $(TEX_KEY)
	cd docs && $(BIB) $(TEX_KEY)

docs$(DIR_SEP)$(TEX_KEY).txt: docs$(DIR_SEP)$(TEX_KEY).pdf
	$(PDFTOTEXT) docs$(DIR_SEP)$(TEX_KEY).pdf docs$(DIR_SEP)$(TEX_KEY).txt

docs$(DIR_SEP)blastersim-out-1.txt: blastersim$(BINEXT)
	$(RUN)blastersim$(BINEXT) > docs$(DIR_SEP)blastersim-out-1.txt

docs$(DIR_SEP)blastersim-out-2.txt: blastersim$(BINEXT) docs$(DIR_SEP)springer-example.nml
	$(RUN)blastersim$(BINEXT) docs$(DIR_SEP)springer-example.nml > docs$(DIR_SEP)blastersim-out-2.txt

# Why move these files here? I don't know how to go up a directory in a cross-platform way in LaTeX. So putting everything in this directory is the easiest approach.
docs$(DIR_SEP)rev.tex: rev.tex
	$(CP) rev.tex docs$(DIR_SEP)rev.tex

docs$(DIR_SEP)test_exact.tex: test_exact.tex
	$(CP) test_exact.tex docs$(DIR_SEP)test_exact.tex

docs$(DIR_SEP)geninput_springer.tex: src$(DIR_SEP)geninput_springer.tex
	$(CP) src$(DIR_SEP)geninput_springer.tex docs$(DIR_SEP)geninput_springer.tex

docs$(DIR_SEP)springer-example.nml: examples$(DIR_SEP)springer-example.nml
	$(CP) examples$(DIR_SEP)springer-example.nml docs$(DIR_SEP)springer-example.nml

docs$(DIR_SEP)defaults.tex: defaults.tex
	$(CP) defaults.tex docs$(DIR_SEP)defaults.tex

# <https://math.nist.gov/~BMiller/LaTeXML/manual/commands/latexml.html>
# Spell checking all HTML files is commented out as it's hard to get aspell to skip code blocks.
# It appears that I can skip by HTML tag: <http://aspell.net/man-html/The-Options.html>
# I am checking the bibliography, however, as aspell won't check that in the TeX files as it's not in a TeX file.
docs$(DIR_SEP)index.html: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bib $(TEX_DEPS)
	cd docs && latexmlc --strict --split --splitnaming=label --dest=index.html $(TEX_KEY).tex
	#$(LOOP_START) docs$(DIR_SEP)*.html $(LOOP_MIDDLE) $(SPELL_HTML) $(LOOP_END)
	$(SPELL_HTML) docs$(DIR_SEP)bib.html
	-$(GREP) "$(BL)Warning$(BR)" docs$(DIR_SEP)$(TEX_KEY).latexml.log

.PHONY: clean_tex
clean_tex:
	$(RM) $(CLEAN_TEX)

# Why not put `$(SPELL_TEX) $(TEX_KEY).tex` in `$(TEX_KEY).bbl`? I get the error "Error: Stdin not a terminal." if I do that.
# TODO: Validate HTML and CSS
# <https://github.com/dspinellis/latex-advice#latex-formatting>: > Users can define their own rules, too, for example to enforce consistency in the spelling of certain phrases. Now it is up to you to create ChkTeX rules to enforce the guidelines in this document. If you do, please share them.
# TODO: Try other BibTeX linters.
# TODO: <https://github.com/sylvainhalle/textidote>
# TODO: /home/ben/git/bmtreport/diction/
.PHONY: lint_tex
lint_tex: docs$(DIR_SEP)$(TEX_KEY).tex docs$(DIR_SEP)$(TEX_KEY).bib docs$(DIR_SEP)index.html docs$(DIR_SEP)$(TEX_KEY).txt docs$(DIR_SEP)$(TEX_KEY).log
	-$(GREP) ?? docs$(DIR_SEP)$(TEX_KEY).txt
	-$(GREP) "$(BL)Underfull \\\hbox$(BR)" docs$(DIR_SEP)$(TEX_KEY).log
	#diction --beginner --suggest docs$(DIR_SEP)$(TEX_KEY).txt
	#style --print-ari 20 docs$(DIR_SEP)$(TEX_KEY).txt
	linkchecker --check-extern docs$(DIR_SEP)index.html
