#######################
# Manual dependencies #
#######################

src$(DIR_SEP)rev.f90: $(ALLSRC)
	$(PYTHON) py$(DIR_SEP)gitrev.py

# jom needs explicit dependencies listed for all files.
src$(DIR_SEP)rev.$(OBJEXT): src$(DIR_SEP)rev.f90

###########################
# Files to manually clean #
###########################

CLEAN_MANUAL = 
