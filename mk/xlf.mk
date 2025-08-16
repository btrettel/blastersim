FFLAGS   = -qlanglvl=2008pure -qwarn64 -qinfo=all -qsuppress=1512-045:1521-004:1516-329 -qsclk=micro
DFLAGS   = -g2 -qcheck=all -qflttrap=inexact:invalid:overflow:underflow:zerodivide
RFLAGS   = -O2 -qreport=hotlist -qlto
AFLAGS   = 
NFLAGS   = -qarch=qauto
SFLAGS   = 
OMPFLAGS = -qsmp=omp

# `-qsuppress` explanations:
# `1512-045`: IBM XL Fortran doesn't seem to support `g0` edit descriptors, so I disable the message for this.
# `1521-004`: Who cares about host association?
# `1516-329`: Seems to have false positives for using variables before setting them.

# `-qsclk=micro`: Needed for test `timeit (absolute)` to pass.

# `-qfree=f90`: free form, default for xlf2008
# `-qsmp=omp`: Will this enable parallelization for `do concurrent`?
# `-qhalt=i`: I would like to include this, but reducing the level to even `-qhalt=e` causes FLT to not compile despite not printing any messages!?!

# Some complaints IBM XL Fortran has that I don't think are valid:
# > "app/genunits.f90", line 23.30: 1518-378 (L) The presence of dummy argument value corresponding to this actual argument might differ from the presence specified by the Fortran 2008 standard.
# Based on my reading of [this IBM forum post](https://community.ibm.com/community/user/power/blogs/archive-user/2013/10/18/using-qxlf2008-to-conform-to-fortran-2008-standard-on-dummy-argument-presence-checking), my usage here is okay.

# TODO: static linking
# <https://www.ibm.com/docs/en/openxl-fortran-aix/17.1.2?topic=programs-dynamic-static-linking>
# > You can use `-b` linker options on the compiler command line to create statically linked object files:

# <https://www.ibm.com/docs/en/xl-fortran-linux/16.1.1?topic=options-qsmp>
