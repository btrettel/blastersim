FFLAGS   = -hlist=l -herror_on_warning -ec
DFLAGS   = -G0 -O0 -hfp0 -Rbcdps
RFLAGS   = -hlist=m -O2 -hfp3
AFLAGS   = 
NFLAGS   = 
SFLAGS   = -hstatic
OMPFLAGS = -homp

# What is `-h thread_do_concurrent`?
# <https://fortran-lang.discourse.group/t/do-concurrent-compiler-flags-to-enable-parallelization/4300/6?u=btrettel>

# `-fsanitize=address`: Seems to have false positives.
# `-M1077`: Disables a warning if OpenMP directives are used but OpenMP is not enabled.

# `-R`
# `b`: > Enables checking of array bounds. Bounds checking is not performed on arrays dimensioned as (1). Enables -Ooverindex.
# `c`: > Enables conformance checking of array operands in array expressions.
# `d`: > Enables a run time check for the `!dir$ collapse` directive and checks the validity of the loop_info count information.
# `p`: > Generates run time code to check the association or allocation status of referenced `POINTER` variables, `ALLOCATABLE` arrays, or assumed-shape arrays.
# `s`: > Enables checking of character substring bounds.
