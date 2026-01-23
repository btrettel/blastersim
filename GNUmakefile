# Summary: Makefile for GNU Make
# Author: Ben Trettel (<http://trettel.us/>)
# Project: [BlasterSim](https://github.com/btrettel/blastersim)
# License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

# Tested with GNU Make on Linux. Works with `make` on Ubuntu.

.POSIX:

# non-POSIX
# <https://innolitics.com/articles/make-delete-on-error/>
.DELETE_ON_ERROR:
MAKEFLAGS = --warn-undefined-variables
# Later, when Make 4.4 is available, add `--shuffle` to `MAKEFLAGS`.

# Detect the operating system so that this can work under both Windows and Linux.
# <https://stackoverflow.com/a/12099167/1124489>
OS ?=
ifeq ($(OS),Windows_NT)
OS_FILENAME = windows
$(error GNU Make building is disabled at the moment on Windows due to weird problems. Try NMAKE or jom.)
# Why use NMAKE or jom on Windows? Some things about GNU Make don't work as well. [Windows directory separators require escaping in a weird way], and [sometimes you need to specify the shell on Windows](https://stackoverflow.com/a/47896799/1124489) as GNU Make will want to use sh.exe. And if you do specify the shell, GNU Make will exit into a new shell.
else
# Despite saying Linux here, I expect this to work on macOS and BSD.
OS_FILENAME = linux
include mk/linux_defaults.mk
endif

# defaults
ifdef F90
FC=$(F90)
else
ifndef FC
FC = gfortran
# fort77 is GNU Make's default, which I'm overriding.
else ifeq ($(FC),fort77)
FC = gfortran
endif
endif
BUILD  = debug
OPENMP = no
STATIC = no

#############
# Compilers #
#############

ifeq ($(FC),gfortran)
include mk$(DIR_SEP)gfortran.mk
else ifeq ($(FC),ifx)
include mk$(DIR_SEP)ifx_$(OS_FILENAME).mk
else ifeq ($(FC),ifort)
include mk$(DIR_SEP)ifort_$(OS_FILENAME).mk
else ifeq ($(FC),nvfortran)
include mk$(DIR_SEP)nvfortran.mk
else ifeq ($(FC),lfortran)
include mk$(DIR_SEP)lfortran.mk
else ifeq ($(FC),crayftn)
include mk$(DIR_SEP)crayftn.mk
else ifeq ($(FC),ftn)
include mk$(DIR_SEP)crayftn.mk
else ifeq ($(FC),xlf2008)
include mk$(DIR_SEP)xlf.mk
else ifeq ($(FC),nagfor)
include mk$(DIR_SEP)nagfor.mk
else ifeq ($(FC),flang-new-19)
include mk$(DIR_SEP)flang.mk
#else
#$(error Invalid FC: $(FC))
endif

ifeq ($(ARCH),native)
AFLAGS = $(NFLAGS)
endif

ifeq ($(BUILD),debug)
FFLAGS += $(DFLAGS) $(AFLAGS)
else ifeq ($(BUILD),release)
FFLAGS += $(RFLAGS) $(AFLAGS)
else
$(error Set BUILD to either debug or release. BUILD=$(BUILD))
endif

ifeq ($(OPENMP),yes)
FFLAGS += $(OMPFLAGS)
else ifeq ($(OPENMP),no)
FFLAGS += 
else
$(error Set OPENMP to either yes or no. OPENMP=$(OPENMP))
endif

ifeq ($(STATIC),yes)
FFLAGS += $(SFLAGS)
else ifeq ($(STATIC),no)
FFLAGS += 
else
$(error Set STATIC to either yes or no. STATIC=$(STATIC))
endif

.f90.$(OBJEXT):
	$(FC) $(OBJFLAGS) $@ $(FFLAGS) $<

include mk$(DIR_SEP)before.mk
include mk$(DIR_SEP)common.mk
include mk$(DIR_SEP)manual.mk
include mk$(DIR_SEP)depends.mk
include mk$(DIR_SEP)latex.mk
include mk$(DIR_SEP)$(OS_FILENAME)_2.mk
